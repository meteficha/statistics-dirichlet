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
    ,empty
    ,Component
    ,fromList
    ,toList
    ,fromDD
    -- * Functions
    ,derive
    ,cost)
    where

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
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
    DM {dmWeights   :: {-# UNPACK #-} !(U.Vector Double)
       ,dmDensities :: {-# UNPACK #-} !(V.Vector DirichletDensity)}
        deriving (Eq)


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

-- | @empty q n x@ is an \"empty\" Dirichlet mixture with @q@
-- components.  Each component has size @n@, weight inversely
-- proportional to its index and all alphas set to @x@.
empty :: Int -> Int -> Double -> DirichletMixture
empty q n x = let dd   = D.empty n x
                  f i  = fromIntegral (q-i) / sum_
                  sum_ = fromIntegral (q*(q+1)`div`2)
              in DM {dmWeights   = U.generate q f
                    ,dmDensities = V.replicate q dd}
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
      qs = U.fromList $ map               fst  components
      as = V.fromList $ map (D.fromList . snd) components

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
       _               -> DM qs as

-- | @toList dm@ is the inverse of @fromList@, constructs a list
-- of components from a Dirichlet mixture.  There are no error
-- conditions and @toList . fromList == id@.
toList :: DirichletMixture -> [Component]
toList (DM qs as) =
    let qs' = U.toList qs
        as' = V.toList $ V.map D.toList as
    in zip qs' as'

-- | Constructs a Dirichlet mixture of one component from a
-- Dirichlet density.
fromDD :: DirichletDensity -> DirichletMixture
fromDD = DM (U.singleton 1) . V.singleton




-- | /Prob(a_j | n, theta)/ Defined in equation (16), "the
-- posterior probability of the /j/-th component of the mixture
-- given the vector of counts /n/".  We return the probabilities
-- for all /j/ in each vector.
--
-- The order of the result is inversed for performance.  In the
-- outer boxed vector there are /j/ elements.  The /i/-th inner
-- unboxed vector contains that probability for each of the
-- training vectors.
--
-- Calculated as per equation (39) using 'logBeta'.
prob_a_n_theta :: TrainingVectors -> DirichletMixture -> V.Vector (U.Vector Double)
prob_a_n_theta ns (DM qs as) =
    let -- Precalculate logBeta of all components
        !logBetaAlphas  = G.unstream $ G.stream $ V.map (logBeta . unDD) as

        -- Calculate the factors for one of the training vectors.
        calc n i q lb_a = let a = unDD (as V.! i)
                          in q * exp (logBeta (U.zipWith (+) n a) - lb_a)
        factors n       = let fs = U.izipWith (calc n) qs logBetaAlphas
                              total = U.sum fs
                          in U.map (/ total) fs
    in transpose (V.length as) $ V.map factors ns


-- | Customized version of @prob_a_n_theta@ used when the weights
-- are being estimated.  Precomputes everything that doesn't
-- depend on the weight.
prob_a_n_theta_w :: TrainingVectors -> V.Vector DirichletDensity
                 -> (U.Vector Double -> V.Vector (U.Vector Double))
prob_a_n_theta_w ns as =
    let -- Precalculate logBeta of all components
        !logBetaAlphas   = G.unstream $ G.stream $ V.map (logBeta . unDD) as

        -- Precalculate the factors for one of the training vectors.
        precalc n i lb_a = let a = unDD (as V.! i)
                           in exp (logBeta (U.zipWith (+) n a) - lb_a)
        !prefactors      = V.map (\n -> U.imap (precalc n) logBetaAlphas) ns
        !as_length       = V.length as

    in \qs ->
        let -- Calculate the final factors.
            calc pfs = let fs = U.zipWith (*) pfs qs
                           total = U.sum fs
                       in U.map (/ total) fs
        in transpose as_length $ V.map calc prefactors


transpose :: Int -> V.Vector (U.Vector Double) -> V.Vector (U.Vector Double)
transpose !as_length !vs =
    V.generate as_length $ \i -> G.unstream $ G.stream $ V.map (U.! i) vs


-- | Cost function for deriving a Dirichlet mixture (equation
-- 18).  This function is minimized by 'derive'.  Calculated
-- using (17) and (54).
cost :: TrainingVectors -> DirichletMixture -> Double
cost ns dm@(DM _ as) =
    let ns_sums    = G.unstream $ G.stream $ V.map U.sum ns
        as_sums    = G.unstream $ G.stream $ V.map (U.sum . unDD) as
    in costWorker (ns, ns_sums) dm as_sums


-- | Worker of 'cost' function that avoids repeating some
-- computations that are done in when reestimating alphas.
costWorker :: (TrainingVectors, U.Vector Double) -> DirichletMixture -> U.Vector Double -> Double
costWorker (!ns, !ns_sums) (DM !qs !as) !as_sums =
    let -- From the equation (54).
        prob_n_a !n !n_sum !a !a_sum !lngamma_a_sum =
            let !s = lngamma (n_sum+1) + lngamma_a_sum - lngamma (n_sum+a_sum)
                f n_i a_i = lngamma (n_i + a_i) - lngamma (n_i + 1) - lngamma a_i
            in exp $ s + U.sum (U.zipWith f n a)

        -- From equation (17).
        prob_n_theta i n =
            let !n_sum = ns_sums U.! i
            in U.sum $ U.zipWith (*) qs $
               U.izipWith (prob_n_a n n_sum . unDD . (as V.!))
                  as_sums lngamma_as_sums
        !lngamma_as_sums = U.map lngamma as_sums
    in negate $ V.sum $ V.imap ((log .) . prob_n_theta) ns

-- | Version of 'cost' function that avoids repeating a lot of
-- computations that are done when reestimating weights.
costWeight :: (TrainingVectors, U.Vector Double) -> V.Vector DirichletDensity
           -> U.Vector Double -> (U.Vector Double -> Double)
costWeight (!ns, !ns_sums) !as !as_sums =
    let -- From the equation (54).
        prob_n_a !n !n_sum !a !a_sum !lngamma_a_sum =
            let !s = lngamma (n_sum+1) + lngamma_a_sum - lngamma (n_sum+a_sum)
                f n_i a_i = lngamma (n_i + a_i) - lngamma (n_i + 1) - lngamma a_i
            in exp $ s + U.sum (U.zipWith f n a)

        -- From equation (17).
        prepare_prob_n_theta i n =
            let !n_sum = ns_sums U.! i
            in {- U.sum $ U.zipWith (*) qs $ -}
               U.izipWith (prob_n_a n n_sum . unDD . (as V.!))
                  as_sums lngamma_as_sums
        !lngamma_as_sums = U.map lngamma as_sums
        !prepared = V.imap prepare_prob_n_theta ns

        -- Final worker function.
        final qs = log . U.sum . U.zipWith (*) qs
    in \(!qs) -> negate $ V.sum $ V.map (final qs) prepared

-- | Derive a Dirichlet mixture using a maximum likelihood method
-- as described by Karplus et al (equation 25).  All training
-- vectors should have the same length, however this is not
-- verified.
derive :: DirichletMixture -> Predicate -> StepSize
         -> TrainingVectors -> Result DirichletMixture
derive (DM initial_qs initial_as)
       (Pred maxIter' minDelta_ deltaSteps' jumpDelta_)
       (Step step) ns
    | V.length ns == 0        = err "empty training data"
    | U.length initial_qs < 1 = err "empty initial weights vector"
    | V.length initial_as < 1 = err "empty initial alphas vector"
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

      -- Compensate the different deltaSteps.
      !minDelta'    = minDelta_  * fromIntegral deltaSteps'
      !jumpDelta'   = jumpDelta_ * fromIntegral deltaSteps'

      -- Sums of each training sequence.
      ns_sums :: U.Vector Double
      !ns_sums = G.unstream $ G.stream $ V.map U.sum ns

      -- Reciprocal of the number of training sequences.
      !recip_m = recip $ fromIntegral $ V.length ns

      -- Functions that work on the alphas only (and not their logs).
      calc_as_sums = G.unstream . G.stream . V.map (U.sum . unDD)

      -- Start training in the zero-th iteration and with
      -- infinite inital cost.
      train =
        let ws      = V.map (U.map log . unDD) as
            as      = initial_as
            as_sums = calc_as_sums as
        in trainAlphas 1 infinity initial_qs ws as as_sums

      trainAlphas !iter !oldCost !qs !ws !as !as_sums =
        -- Calculate Prob(a | n, theta)
        let !probs_a_n   = prob_a_n_theta ns (DM qs as)

        -- Calculate all S_j's.
            !sjs         = G.unstream $ G.stream $ V.map U.sum probs_a_n

        -- Reestimate alpha's.
            calc_ws !j =
              -- Everything that doesn't depend on i, just on j.
              let a_sum        = as_sums U.! j
                  psi_a_sum    = psi a_sum
                  probs        = probs_a_n V.! j
                  sum_prob_psi = U.sum $ U.zipWith (*) probs $
                                 U.map (psi . (+) a_sum) ns_sums
              -----
              in \(!i) !w_i !a_i ->
                let !s1 = (sjs U.! j) * (psi_a_sum - psi a_i)
                    !s2 = V.sum $ V.map f ns
                    f n = (probs U.! i) * psi (n U.! i + a_i)
                in w_i + step * a_i * (s1 + s2 - sum_prob_psi)
            !ws' = V.izipWith (\j w (DD a) -> U.izipWith (calc_ws j) w a) ws as
            !as' = V.map (DD . U.map exp) ws'

        -- Recalculate constants.
            !as_sums' = calc_as_sums as'
            !calcCost = iter `mod` deltaSteps' == 0
            !cost'    = if calcCost then newCost else oldCost
             where newCost = costWorker (ns, ns_sums) (DM qs as') as_sums'
            !delta    = abs (cost' - oldCost)

        -- Verify convergence.  Even with MaxIter we only stop
        -- iterating if the delta was calculated.  Otherwise we
        -- wouldn't be able to tell the caller why the delta was
        -- still big when we reached the limit.
        in case (calcCost, delta <= minDelta', delta <= jumpDelta', iter >= maxIter') of
             (True,True,_,_) -> Result Delta   iter delta cost' $ DM qs as'
             (True,_,True,_) -> trainWeights (iter+1) cost' qs ws' as' as_sums'
             (True,_,_,True) -> Result MaxIter iter delta cost' $ DM qs as'
             _               -> trainAlphas  (iter+1) cost' qs ws' as' as_sums'

      trainWeights !oldIter !veryOldCost !oldQs !ws !as !as_sums =
        -- Prepare invariant parts.
        let !probs_a_n_mk = prob_a_n_theta_w ns as
            !cost_mk      = costWeight (ns, ns_sums) as as_sums
        in ($ oldQs) . ($ veryOldCost) . ($ oldIter) . fix $ \again !iter !oldCost !qs ->
          -- Reestimate weight's.
          let !probs_a_n = probs_a_n_mk qs
              qs' = G.unstream $ G.stream $ V.map ((*) recip_m . U.sum) probs_a_n

          -- Recalculate constants.
              !calcCost = iter `mod` deltaSteps' == 0
              !cost'    = if calcCost then cost_mk qs' else oldCost
              !delta    = abs (cost' - oldCost)

        -- Verify convergence.  We never stop the process here.
        in case (calcCost, delta <= jumpDelta') of
             (True,True) -> trainAlphas (iter+1) cost' qs' ws as as_sums
             _           -> again (iter+1) cost' qs'



