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
    ,dmParameters
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
import qualified Data.Vector.Unboxed as U

import Control.Parallel.Strategies (NFData(..), rwhnf)
import Control.Monad.ST
import Data.Function (fix)
import Numeric.GSL.Special.Gamma (lngamma)
import Numeric.GSL.Special.Psi (psi)

import qualified Numeric.Optimization.Algorithms.HagerZhang05 as CG

import qualified Math.Statistics.Dirichlet.Density as D
import qualified Math.Statistics.Dirichlet.Matrix as M
import Math.Statistics.Dirichlet.Density (DirichletDensity(..))
import Math.Statistics.Dirichlet.Matrix (Matrix (..))
import Math.Statistics.Dirichlet.Options
import Math.Statistics.Dirichlet.Util



-- | A Dirichlet mixture.
data DirichletMixture =
    DM {dmWeights    :: {-# UNPACK #-} !(U.Vector Double)
        -- ^ Weights of each density.
       ,dmDensities  :: {-# UNPACK #-} !M.Matrix
        -- ^ Values of all parameters of all densities.  This
        -- matrix has @length dmWeights@ rows.
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

-- | Number of parameters each component has.
dmParameters :: DirichletMixture -> Int
dmParameters = mCols . dmDensities

-- | Separated list of densities.
dmDensitiesL :: DirichletMixture -> [DirichletDensity]
dmDensitiesL (DM _ as) = map DD $ V.toList $ M.rows as

-- | @dm !!! i@ is the @i@-th density.  No bounding checks are
-- made.
(!!!) :: DirichletMixture -> Int -> U.Vector Double
(DM _ as) !!! i = as M.!!! i
{-# INLINE (!!!) #-}





dmap :: (U.Vector Double -> Double) -> DirichletMixture -> U.Vector Double
dmap f = M.rowmap f . dmDensities



-- | @empty q n x@ is an \"empty\" Dirichlet mixture with @q@
-- components and @n@ parameters.  Each component has size @n@,
-- weight inversely proportional to its index and all alphas set
-- to @x@.
empty :: Int -> Int -> Double -> DirichletMixture
empty q n x = let (DD d) = D.empty n x
                  f i    = fromIntegral (q-i) / sum_
                  sum_   = fromIntegral (q*(q+1)`div`2)
              in DM {dmWeights    = U.generate q f
                    ,dmDensities  = M.replicateRows q d}
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
      qs =         U.fromList $       map fst components
      as = M q n $ U.fromList $ concatMap snd components

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
toList dm =
    let qs' = U.toList $ dmWeights dm
        as' = map (U.toList . unDD) (dmDensitiesL dm)
    in zip qs' as'

-- | Constructs a Dirichlet mixture of one component from a
-- Dirichlet density.
fromDD :: DirichletDensity -> DirichletMixture
fromDD (DD d) = DM (U.singleton 1) (M.replicateRows 1 d)











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
prob_a_n_theta :: TrainingVectors -> DirichletMixture -> Matrix
prob_a_n_theta ns dm@(DM qs _) =
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
    in M.fromVectorT $ V.map factors ns


-- | Customized version of @prob_a_n_theta@ used when the weights
-- are being estimated.  Precomputes everything that doesn't
-- depend on the weight.
prob_a_n_theta_weights :: TrainingVectors -> Matrix
                       -> (U.Vector Double -> Matrix)
prob_a_n_theta_weights ns as =
    let -- Precalculate logBeta of all components
        !logBetaAlphas   = M.rowmap logBeta as

        -- Precalculate the factors for one of the training vectors.
        precalc n i lb_a = let !a = as M.!!! i
                           in logBeta (U.zipWith (+) n a) - lb_a
        norm fs          = let !c = U.maximum fs
                           in U.map (exp . subtract c) fs
        !prefactors      = V.map (norm . flip U.imap logBetaAlphas . precalc) ns

    in \qs ->
        let -- Calculate the final factors.
            calc pfs = let fs = U.zipWith (*) pfs qs
                           total = U.sum fs
                       in U.map (/ total) fs
        in M.fromVectorT $ V.map calc prefactors












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
cost_worker (!ns, !ns_sums) dm@(DM !qs _) !as_sums =
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
cost_weight :: (TrainingVectors, U.Vector Double) -> Matrix
            -> U.Vector Double -> (U.Vector Double -> Double)
cost_weight (!ns, !ns_sums) !as !as_sums =
    let -- From the equation (54).
        prob_n_a !n !n_sum !a !a_sum !lngamma_a_sum =
            let !s = lngamma (n_sum+1) + lngamma_a_sum - lngamma (n_sum+a_sum)
                f n_i a_i = lngamma (n_i + a_i) - lngamma (n_i + 1) - lngamma a_i
            in exp $ s + U.sum (U.zipWith f n a)

        -- From equation (17).
        prepare_prob_n_theta i n =
            let !n_sum = ns_sums U.! i
            in {- U.sum $ U.zipWith (*) qs $ -}
               U.izipWith (prob_n_a n n_sum . (as M.!!!))
                  as_sums lngamma_as_sums
        !lngamma_as_sums = U.map lngamma as_sums
        !prepared = V.imap prepare_prob_n_theta ns

        -- Final worker function.
        final qs = log . U.sum . U.zipWith (*) qs
    in \(!qs) -> negate $ V.sum $ V.map (final qs) prepared







-- | Derivative of the cost function with respect @w_{i,j}@,
-- defined by Equation (22).  The result is given in the same
-- size and order as the 'dmDensitites' vector.
del_cost_w :: TrainingVectors -> DirichletMixture -> Matrix
del_cost_w ns dm =
    let ns_sums = G.unstream $ G.stream $ V.map U.sum ns
        as_sums = dmap U.sum dm
        tns     = M.fromVectorT ns
    in del_cost_w_worker (ns, tns, ns_sums) dm as_sums


-- | Worker function of 'del_cost_w'.
del_cost_w_worker :: (TrainingVectors, Matrix, U.Vector Double)
                  -> DirichletMixture -> U.Vector Double -> Matrix
del_cost_w_worker (!ns, !tns, !ns_sums) dm !as_sums =
    let -- Calculate Prob(a | n, theta)
        !probs_a_n   = prob_a_n_theta ns dm

        -- Calculate all S_j's.
        !sjs         = M.rowmap U.sum probs_a_n

        -- @calc j _ i _ _@ calculates the derivative of the
        -- cost function with respect to @w_{i,j}@.  The other
        -- arguments come from vector that we @zipWith@ below.
        calc j probs =
          -- Everything that doesn't depend on i, just on j.
          let !a_sum        = as_sums U.! j
              !psi_a_sum    = psi a_sum
              !sum_prob_psi = U.sum $ U.zipWith (*) probs $
                              U.map (psi . (+) a_sum) ns_sums
          -----
          in \i a_i ->
            let !s1 = (sjs U.! j) * (psi_a_sum - psi a_i)
                !s2 = U.sum $ U.zipWith (\p_i n_i -> p_i * psi (n_i + a_i)) probs (tns M.!!! i)
            in - a_i * (s1 + s2 - sum_prob_psi)

    in M.fromVector $ V.imap (\j p_j -> let !f = calc j p_j
                                        in U.imap f (dm !!! j))
                             (M.rows probs_a_n)





-- | Derive a Dirichlet mixture using a maximum likelihood method
-- as described by Karplus et al (equation 25) using CG_DESCENT
-- method by Hager and Zhang (see
-- "Numeric.Optimization.Algorithms.HagerZhang05").  All training
-- vectors should have the same length, however this is not
-- verified.
derive :: DirichletMixture -> Predicate -> StepSize
         -> TrainingVectors -> Result DirichletMixture
derive (DM initial_qs initial_as)
       (Pred maxIter minDelta _deltaSteps maxWeightIter jumpDelta)
       (Step step) ns
    | V.length ns == 0          = err "empty training data"
    | U.length initial_qs < 1   = err "empty initial weights vector"
    | M.size initial_as < (1,1) = err "empty initial alphas vector"
    | maxIter < 1               = err "non-positive maxIter"
    | minDelta < 0              = err "negative minDelta"
    | jumpDelta < 0             = err "negative jumpDelta"
    | jumpDelta < minDelta      = err "minDelta greater than jumpDelta"
    | step <= 0                 = err "non-positive step"
    | step >= 1                 = err "step greater than one"
    | otherwise                 = runST train
    where
      err = error . ("Dirichlet.derive: " ++)
      singleDensity = U.length initial_qs == 1

      -- Sums of each training sequence.
      ns_sums :: U.Vector Double
      !ns_sums = G.unstream $ G.stream $ V.map U.sum ns

      -- Transposed training sequences.
      !tns = M.fromVectorT ns

      -- Reciprocal of the number of training sequences.
      !recip_m = recip $ fromIntegral $ V.length ns

      -- Functions that work on the alphas only (and not their logs).
      calc_as_sums = M.rowmap U.sum

      -- Parameters used by CG_DESCENT.
      parameters = CG.defaultParameters
                     {CG.printFinal    = True
                     ,CG.printParams   = True
                     ,CG.verbose       = CG.VeryVerbose
                     ,CG.maxItersFac   = max 1 $ fromIntegral maxIter / 20
                     ,CG.estimateError = CG.RelativeEpsilon (1e-6 * s)}
        where (w,h) = M.size initial_as
              s = fromIntegral (w * h * V.length ns)

      -- Transform a U.Vector from/to a M.Matrix in the case that
      -- the matrix has the same shape as initial_as (i.e. all
      -- as's and ws's).
      fromMatrix = M.mData
      toMatrix v = initial_as {M.mData = v}

      -- Create specialized functions that are optimized by
      -- CG_DESCENT.  They depend only on @qs@, the weights.
      createFunctions !qs =
        let calc f = \ws -> let !as      = M.map exp (toMatrix ws)
                                !as_sums = calc_as_sums as
                                dm       = DM qs as
                            in f dm as_sums
            grad_worker = ((fromMatrix .) .) . del_cost_w_worker
            func = CG.VFunction $ calc $ cost_worker (ns, ns_sums)
            grad = CG.VGradient $ calc $ grad_worker (ns, tns, ns_sums)
            comb = CG.VCombined $ calc $ \dm as_sums ->
                     (cost_worker (ns, ns_sums) dm as_sums
                     ,grad_worker (ns, tns, ns_sums) dm as_sums)
        in (func, grad, comb)

      -- Start training in the zero-th iteration and with
      -- infinite inital cost.
      train = trainAlphas 0 infinity initial_qs (M.map log initial_as)

      trainAlphas !iter !oldCost !qs !ws = {-# SCC "trainAlphas" #-} do
        -- Optimize using CG_DESCENT
        let (func, grad, comb) = createFunctions qs
            opt = CG.optimize parameters minDelta (fromMatrix ws)
                              func grad (Just comb)

        (!pre_ws', result, stats) <- unsafeIOToST opt
        let !ws' = toMatrix (G.unstream $ G.stream pre_ws')

        -- Recalculate everything.
            !as'      = M.map exp ws'
            !as_sums' = calc_as_sums as'
            !iter'    = iter + fromIntegral (CG.totalIters stats)
            !cost'    = CG.finalValue stats
            !delta    = abs (cost' - oldCost)

        -- Verify convergence.  Even with MaxIter we only stop
        -- iterating if the delta was calculated.  Otherwise we
        -- wouldn't be able to tell the caller why the delta was
        -- still big when we reached the limit.
        case (decide result
             ,delta <= minDelta
             ,iter' >= maxIter
             ,singleDensity) of
            (Error,_,_,_)  -> error $ "Mixture.derive: CG_DESCENT returned " ++ show result
            (Stop r,_,_,_) -> return $ Result r       iter' delta cost' $ DM qs as'
            (_,True,_,_)   -> return $ Result Delta   iter' delta cost' $ DM qs as'
            (_,_,True,_)   -> return $ Result MaxIter iter' delta cost' $ DM qs as'
            (_,_,_,True)   -> return $ Result Delta   iter' delta cost' $ DM qs as'
            (GoOn,_,_,_)   -> trainWeights iter' cost' qs ws' as' as_sums'

      trainWeights !oldIter !veryOldCost !oldQs !ws !as !as_sums =
        {-# SCC "trainWeights" #-}
        -- Prepare invariant parts.
        let !probs_a_n_mk = prob_a_n_theta_weights ns as
            !cost_mk      = cost_weight (ns, ns_sums) as as_sums
        in ($ oldQs) . ($ veryOldCost) . ($ maxWeightIter) . fix $
               \again !itersLeft !oldCost !qs ->
          -- Reestimate weight's.
          let !probs_a_n = probs_a_n_mk qs
              qs' = M.rowmap ((*) recip_m . U.sum) probs_a_n

          -- Recalculate constants.
              !cost'    = cost_mk qs'
              !delta    = abs (cost' - oldCost)

        -- Verify convergence.  We never stop the process here.
        in case (delta <= jumpDelta, itersLeft <= 0) of
             (False,False) -> again (itersLeft-1) cost' qs'
             _             -> trainAlphas oldIter cost' qs' ws


-- | Decide what we should do depending on the result of the
-- CG_DESCENT routine.
decide :: CG.Result -> Decision
decide CG.ToleranceStatisfied      = GoOn
decide CG.FunctionChange           = GoOn
decide CG.MaxTotalIter             = GoOn
decide CG.MaxSecantIter            = GoOn

decide CG.FunctionValueNaN         = Stop Delta
decide CG.LineSearchFailsInitial   = Stop Delta
decide CG.LineSearchFailsBisection = Stop Delta
decide CG.LineSearchFailsUpdate    = Stop Delta

decide CG.StartFunctionValueNaN    = Error
decide CG.NegativeSlope            = Error
decide CG.NotDescent               = Error
decide CG.DebugTol                 = Error

data Decision = GoOn | Stop Reason | Error