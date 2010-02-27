---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet.Mixture
-- Copyright   : (c) 2009 Felipe Lessa
-- License     : GPL
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
--------------------------------------------------------------------------

module Math.Statistics.Dirichlet.Mixture
    (-- * Data types
     DirichletMixture(..)
    ,dmComponents
    ,dmDensitiesL
    ,(!!!)
    ,empty
    ,Component
    ,fromList
    ,toList
    ,fromDD
    -- * Functions
    ,derive
    ,cost
    ,del_cost_w)
    where

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Fusion.Stream as S
import qualified Data.Vector.Unboxed as U

import Control.Parallel.Strategies (NFData(..), rwhnf)
import Data.Function (fix)
import Numeric.GSL.Special.Gamma (lngamma)
import Numeric.GSL.Special.Psi (psi)

import qualified Math.Statistics.Dirichlet.Density as D
import Math.Statistics.Dirichlet.Density (DirichletDensity(..))
import Math.Statistics.Dirichlet.Options
import Math.Statistics.Dirichlet.Util



-- | A Dirichlet mixture.
data DirichletMixture =
    DM {dmParameters :: {-# UNPACK #-} !Int
        -- ^ Number of parameters each density has.
       ,dmWeights    :: {-# UNPACK #-} !(U.Vector Double)
        -- ^ Weights of each density.
       ,dmDensities  :: {-# UNPACK #-} !(U.Vector Double)
        -- ^ Values of all parameters of all densities.  This
        -- vector has @dmParameters * length dmWeights@ values.
       } deriving (Eq)

instance Show DirichletMixture where
    showsPrec prec dm =
      showParen (prec > 10) $
      showString "fromList " .
      showsPrec 11 (toList dm)

instance Read DirichletMixture where
    readsPrec p ('(':xs) = let (ys,')':zs) = break (== ')') xs
                           in map (\(x,s) -> (x,s++zs)) $
                              readsPrec p ys
    readsPrec p xs = let [("fromList",list)] = lex xs
                     in map (\(x,s) -> (fromList x,s)) $
                        readsPrec p list

instance NFData DirichletMixture where
    rnf = rwhnf


-- | Number of components in a dirichlet mixture.
dmComponents :: DirichletMixture -> Int
dmComponents = U.length . dmWeights

-- | Separated list of densities.
dmDensitiesL :: DirichletMixture -> [DirichletDensity]
dmDensitiesL (DM n _ as) =
    let l = U.length as
    in [DD (U.unsafeSlice i n as) | i <- [0,n..l-1]]

-- | @dm !!! i@ is the @i@-th density.  No bounding checks are
-- made.
(!!!) :: DirichletMixture -> Int -> U.Vector Double
(DM n _ as) !!! i = get n as i
{-# INLINE (!!!) #-}





get :: Int -> U.Vector Double -> Int -> U.Vector Double
get n as i = U.unsafeSlice (i*n) n as
{-# INLINE get #-}

dmap :: (U.Vector Double -> Double) -> DirichletMixture -> U.Vector Double
dmap f dm = dmap' f (dmComponents dm) (dmParameters dm) (dmDensities dm)

dmap' :: (U.Vector Double -> Double) -> Int -> Int -> U.Vector Double -> U.Vector Double
dmap' f components parameters as = U.generate components (f . get parameters as)






-- | @empty q n x@ is an \"empty\" Dirichlet mixture with @q@
-- components and @n@ parameters.  Each component has size @n@,
-- weight inversely proportional to its index and all alphas set
-- to @x@.
empty :: Int -> Int -> Double -> DirichletMixture
empty q n x = let (DD d) = D.empty n x
                  f i    = fromIntegral (q-i) / sum_
                  sum_   = fromIntegral (q*(q+1)`div`2)
              in DM {dmParameters = n
                    ,dmWeights    = U.generate q f
                    ,dmDensities  = U.generate (q*n) ((d U.!) . (`mod` n))}
{-# INLINE empty #-}


-- | A list representation of a component of a Dirichlet mixture.
-- Used by 'fromList' and 'toList' only.
type Component = (Double, [Double])

-- | @fromList xs@ constructs a Dirichlet mixture from a
-- non-empty list of components.  Each component has a weight and
-- a list of alpha values.  The weights sum to 1, all lists must
-- have the same number of values and every number must be
-- non-negative.  All of these preconditions are verified for
-- clear mistakes.
fromList :: [Component] -> DirichletMixture
fromList components =
  let -- Vectors
      qs = U.fromList $       map fst components
      as = U.fromList $ concatMap snd components

      -- Properties of the mixture
      q  = length components
      n  = length (snd $ head components)

      -- Checks
      c0 = q >= 1
      c1 = abs (U.sum qs - 1) < 1e-2 -- we're quite permissive here
      c2 = U.all (>= 0) qs
      c3 = all ((== n) . length . snd) components
      c4 = all (all (>= 0)      . snd) components
      e  = error . ("Dirichlet.Mixture.fromList: " ++)
  in case (c0, c1, c2, c3, c4) of
       (False,_,_,_,_) -> e "there must be at least one component"
       (_,False,_,_,_) -> e "the sum of the weights must be one"
       (_,_,False,_,_) -> e "all weights must be greater than or equal to zero"
       (_,_,_,False,_) -> e "every component must have the same size"
       (_,_,_,_,False) -> e "all alphas must be greater than or equal to zero"
       _               -> DM n qs as

-- | @toList dm@ is the inverse of @fromList@, constructs a list
-- of components from a Dirichlet mixture.  There are no error
-- conditions and @toList . fromList == id@.
toList :: DirichletMixture -> [Component]
toList dm@(DM _ qs _) =
    let qs' = U.toList qs
        as' = map (U.toList . unDD) $ dmDensitiesL dm
    in zip qs' as'

-- | Constructs a Dirichlet mixture of one component from a
-- Dirichlet density.
fromDD :: DirichletDensity -> DirichletMixture
fromDD (DD d) = DM (U.length d) (U.singleton 1) d











-- | /Prob(a_j | n, theta)/ Defined in equation (16), "the
-- posterior probability of the j-th component of the mixture
-- given the vector of counts n".  We return the probabilities
-- for all /j/ in each vector.
--
-- The order of the result is inversed for performance.  In the
-- outer boxed vector there are /j/ elements.  The /i/-th inner
-- unboxed vector contains that probability for each of the
-- training vectors.
--
-- Calculated as per equation (39) using 'logBeta'.  If we take
-- the numerator of the right hand side of equation (39) as /Y_j/
-- and the left hand side as /P_j/, then /P_j/ is proportional to
-- /Y_j/ normalized to sum to 1.  We may have problems if /P_j/
-- is too large or too small.  Using the suggestion from the
-- paper, we may multiply all /P_j/ by a constant before
-- normalizing everything.  We calculate /P_j/ using a logarithm,
-- so that means we may freely add or subtract a constant from
-- the logarithm before appling the exponential function.  This
-- is really essencial.
prob_a_n_theta :: TrainingVectors -> DirichletMixture -> V.Vector (U.Vector Double)
prob_a_n_theta ns dm@(DM _ qs _) =
    let -- Precalculate logBeta of all components
        !logBetaAlphas = dmap logBeta dm

        -- Calculate the factors for one of the training vectors.
        calc n i lb_a  = let !a = dm !!! i
                         in logBeta (U.zipWith (+) n a) - lb_a
        factors n      = let fs  = U.imap (calc n) logBetaAlphas
                             !c  = U.maximum fs  -- see the note above
                             fs' = U.zipWith (\q f -> q * exp (f - c)) qs fs
                             !total = U.sum fs'
                         in U.map (/ total) fs'
    in transpose (dmComponents dm) $ V.map factors ns


-- | Customized version of @prob_a_n_theta@ used when the weights
-- are being estimated.  Precomputes everything that doesn't
-- depend on the weight.
prob_a_n_theta_weights :: TrainingVectors -> Int -> Int -> U.Vector Double
                       -> (U.Vector Double -> V.Vector (U.Vector Double))
prob_a_n_theta_weights ns components parameters as =
    let -- Precalculate logBeta of all components
        !logBetaAlphas   = dmap' logBeta components parameters as

        -- Precalculate the factors for one of the training vectors.
        precalc n i lb_a = let !a = get parameters as i
                           in logBeta (U.zipWith (+) n a) - lb_a
        norm fs          = let !c = U.maximum fs
                           in U.map (exp . subtract c) fs
        !prefactors      = V.map (\n -> norm $ U.imap (precalc n) logBetaAlphas) ns

    in \qs ->
        let -- Calculate the final factors.
            calc pfs = let fs = U.zipWith (*) pfs qs
                           total = U.sum fs
                       in U.map (/ total) fs
        in transpose components $ V.map calc prefactors


transpose :: Int -> V.Vector (U.Vector Double) -> V.Vector (U.Vector Double)
transpose !as_length !vs = -- as_length should be equal to U.length (V.head vs)
    V.generate as_length $ \i -> G.unstream $ G.stream $ V.map (U.! i) vs












-- | Cost function for deriving a Dirichlet mixture (equation
-- 18).  This function is minimized by 'derive'.  Calculated
-- using (17) and (54).
cost :: TrainingVectors -> DirichletMixture -> Double
cost ns dm =
    let ns_sums = G.unstream $ G.stream $ V.map U.sum ns
        as_sums = dmap U.sum dm
    in cost_worker (ns, ns_sums) dm as_sums


-- | Worker of 'cost' function that avoids repeating some
-- computations that are done when reestimating alphas.
cost_worker :: (TrainingVectors, U.Vector Double) -> DirichletMixture
            -> U.Vector Double -> Double
cost_worker (!ns, !ns_sums) dm@(DM _ !qs _) !as_sums =
    let -- From the equation (54).
        prob_n_a !n !n_sum !a !a_sum !lngamma_a_sum =
            let !s = lngamma (n_sum+1) + lngamma_a_sum - lngamma (n_sum+a_sum)
                f n_i a_i = lngamma (n_i + a_i) - lngamma (n_i + 1) - lngamma a_i
            in exp $ s + U.sum (U.zipWith f n a)

        -- From equation (17).
        prob_n_theta i n =
            let !n_sum = ns_sums U.! i
            in U.sum $ U.zipWith (*) qs $
               U.izipWith (prob_n_a n n_sum . (dm !!!))
                  as_sums lngamma_as_sums
        !lngamma_as_sums = U.map lngamma as_sums
    in negate $ V.sum $ V.imap ((log .) . prob_n_theta) ns

-- | Version of 'cost' function that avoids repeating a lot of
-- computations that are done when reestimating weights.
cost_weight :: (TrainingVectors, U.Vector Double) -> Int -> U.Vector Double
            -> U.Vector Double -> (U.Vector Double -> Double)
cost_weight (!ns, !ns_sums) !parameters !as !as_sums =
    let -- From the equation (54).
        prob_n_a !n !n_sum !a !a_sum !lngamma_a_sum =
            let !s = lngamma (n_sum+1) + lngamma_a_sum - lngamma (n_sum+a_sum)
                f n_i a_i = lngamma (n_i + a_i) - lngamma (n_i + 1) - lngamma a_i
            in exp $ s + U.sum (U.zipWith f n a)

        -- From equation (17).
        prepare_prob_n_theta i n =
            let !n_sum = ns_sums U.! i
            in {- U.sum $ U.zipWith (*) qs $ -}
               U.izipWith (prob_n_a n n_sum . get parameters as)
                  as_sums lngamma_as_sums
        !lngamma_as_sums = U.map lngamma as_sums
        !prepared = V.imap prepare_prob_n_theta ns

        -- Final worker function.
        final qs = log . U.sum . U.zipWith (*) qs
    in \(!qs) -> negate $ V.sum $ V.map (final qs) prepared







-- | Derivative of the cost function with respect @w_{i,j}@,
-- defined by Equation (22).  The result is given in the same
-- size and order as the 'dmDensitites' vector.
del_cost_w :: TrainingVectors -> DirichletMixture -> U.Vector Double
del_cost_w ns dm =
    let ns_sums = G.unstream $ G.stream $ V.map U.sum ns
        as_sums = dmap U.sum dm
        tns     = transpose (dmComponents dm) ns
    in del_cost_w_worker (ns, tns, ns_sums) dm as_sums


-- | Worker function of 'del_cost_w'.
del_cost_w_worker :: (TrainingVectors, V.Vector (U.Vector Double), U.Vector Double)
                  -> DirichletMixture -> U.Vector Double -> U.Vector Double
del_cost_w_worker (!ns, !tns, !ns_sums) dm !as_sums =
    let -- Calculate Prob(a | n, theta)
        !probs_a_n   = prob_a_n_theta ns dm

        -- Calculate all S_j's.
        !sjs         = G.unstream $ G.stream $ V.map U.sum probs_a_n

        -- @calc j _ _ i _ _@ calculates the derivative of the
        -- cost function with respect to @w_{i,j}@.  The other
        -- arguments come from vector that we @zipWith@ below.
        calc j probs tn_j =
          -- Everything that doesn't depend on i, just on j.
          let !a_sum        = as_sums U.! j
              !psi_a_sum    = psi a_sum
              !sum_prob_psi = U.sum $ U.zipWith (*) probs $
                              U.map (psi . (+) a_sum) ns_sums
          -----
          in \p_i a_i ->
            let !s1 = (sjs U.! j) * (psi_a_sum - psi a_i)
                !s2 = U.sum $ U.map (\n_i -> p_i * psi (n_i + a_i)) tn_j
            in a_i * (s1 + s2 - sum_prob_psi)

    in G.unstream $ S.concatMap G.stream $ G.stream $
       V.izipWith (\j p_j tn_j -> let !f = calc j p_j tn_j
                                  in U.zipWith f p_j (dm !!! j))
                  probs_a_n tns












-- | Derive a Dirichlet mixture using a maximum likelihood method
-- as described by Karplus et al (equation 25).  All training
-- vectors should have the same length, however this is not
-- verified.
derive :: DirichletMixture -> Predicate -> StepSize
         -> TrainingVectors -> Result DirichletMixture
derive idm@(DM parameters initial_qs initial_as)
       (Pred maxIter' minDelta_ deltaSteps' maxWeightIter' jumpDelta_)
       (Step step) ns
    | V.length ns == 0        = err "empty training data"
    | U.length initial_qs < 1 = err "empty initial weights vector"
    | U.length initial_as < 1 = err "empty initial alphas vector"
    | maxIter' < 1            = err "non-positive maxIter"
    | minDelta_ < 0           = err "negative minDelta"
    | jumpDelta_ < 0          = err "negative jumpDelta"
    | jumpDelta_ < minDelta_  = err "minDelta greater than jumpDelta"
    | deltaSteps' < 1         = err "non-positive deltaSteps"
    | step <= 0               = err "non-positive step"
    | step >= 1               = err "step greater than one"
    | otherwise               = train
    where
      err = error . ("Dirichlet.derive: " ++)
      !components = dmComponents idm

      -- Compensate the different deltaSteps.
      !minDelta'    = minDelta_  * fromIntegral deltaSteps'
      !jumpDelta'   = jumpDelta_ * fromIntegral deltaSteps'

      -- Sums of each training sequence.
      ns_sums :: U.Vector Double
      !ns_sums = G.unstream $ G.stream $ V.map U.sum ns

      -- Transposed training sequences.
      !tns = transpose components ns

      -- Reciprocal of the number of training sequences.
      !recip_m = recip $ fromIntegral $ V.length ns

      -- Functions that work on the alphas only (and not their logs).
      calc_as_sums = dmap U.sum . DM parameters initial_qs

      -- Start training in the zero-th iteration and with
      -- infinite inital cost.
      train =
        let ws      = U.map log as
            as      = initial_as
            as_sums = calc_as_sums as
        in trainAlphas 1 infinity initial_qs ws as as_sums

      trainAlphas !iter !oldCost !qs !ws !as !as_sums =
        {-# SCC "trainAlphas" #-}
        -- Calculate derivative, follow in steepest descent.
        let derivative = del_cost_w_worker (ns, tns, ns_sums) (DM parameters qs as) as_sums
            !ws' = U.zipWith ((+) . (step *)) derivative ws
            !as' = U.map exp ws'

        -- Recalculate constants.
            !as_sums' = calc_as_sums as'
            !calcCost = iter `mod` deltaSteps' == 0
            !cost'    = if calcCost then newCost else oldCost
             where newCost = cost_worker (ns, ns_sums) (DM parameters qs as') as_sums'
            !delta    = abs (cost' - oldCost)

        -- Verify convergence.  Even with MaxIter we only stop
        -- iterating if the delta was calculated.  Otherwise we
        -- wouldn't be able to tell the caller why the delta was
        -- still big when we reached the limit.
        in case (calcCost, delta <= minDelta', iter >= maxIter', delta <= jumpDelta') of
             (True,True,_,_) -> Result Delta   iter delta cost' $ DM parameters qs as'
             (True,_,True,_) -> Result MaxIter iter delta cost' $ DM parameters qs as'
             (True,_,_,True) -> trainWeights (iter+1) cost' qs ws' as' as_sums'
             _               -> trainAlphas  (iter+1) cost' qs ws' as' as_sums'

      trainWeights !oldIter !veryOldCost !oldQs !ws !as !as_sums =
        {-# SCC "trainWeights" #-}
        -- Prepare invariant parts.
        let !probs_a_n_mk = prob_a_n_theta_weights ns components parameters as
            !cost_mk      = cost_weight (ns, ns_sums) parameters as as_sums
        in ($ oldQs) . ($ veryOldCost) . ($ maxWeightIter') . fix $
               \again !itersLeft !oldCost !qs ->
          -- Reestimate weight's.
          let !probs_a_n = probs_a_n_mk qs
              qs' = G.unstream $ G.stream $ V.map ((*) recip_m . U.sum) probs_a_n

          -- Recalculate constants.
              !calcCost = itersLeft `mod` deltaSteps' == 0
              !cost'    = if calcCost then cost_mk qs' else oldCost
              !delta    = abs (cost' - oldCost)

        -- Verify convergence.  We never stop the process here.
        in case (calcCost && delta <= jumpDelta', itersLeft <= 0) of
             (False,False) -> again (itersLeft-1) cost' qs'
             _             -> trainAlphas oldIter cost' qs' ws as as_sums



