---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet
-- Copyright   : (c) 2009 Felipe Lessa
-- License     : GPL
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
--------------------------------------------------------------------------

module Math.Statistics.Dirichlet
    (TrainingVector
    ,StepSize(..)
    ,Delta
    ,Predicate(..)
    ,Reason(..)
    ,Result(..)

    ,DirichletDensity(..)
    ,emptyDD
    ,deriveDD
    ,costDD

    ,logBeta)
    where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

import Control.Parallel.Strategies (NFData(..), rwhnf)
import Numeric.GSL.Special.Gamma (lngamma, lnbeta)
import Numeric.GSL.Special.Psi (psi)



-- | Logarithm of the beta function applied to a vector.
logBeta :: U.Vector Double -> Double
logBeta xs | U.length xs == 2 = lnbeta (U.head xs) (U.last xs)
           | otherwise        = U.sum (U.map lngamma xs) - lngamma (U.sum xs)
{-# INLINE logBeta #-}

-- | A vector used for deriving the parameters of a Dirichlet
--   density or mixture.
type TrainingVector = U.Vector Double

-- | Usually denoted by lowercase greek letter eta (η), size of
--   each step in the gradient. Should be greater than zero and
--   much less than one.
newtype StepSize = Step Double

-- | Maximum difference between costs to consider that the
--   process converged.
type Delta = Double

-- | Predicate specifying when the training should be over.
data Predicate = Pred
    {maxIter  :: !Int    -- ^ Maximum number of iterations.
    ,minDelta :: !Delta  -- ^ Minimum delta to continue iterating.
    }
                 deriving (Eq, Read, Show)

-- | Reason why the derivation was over.
data Reason = Delta | MaxIter
              deriving (Eq, Read, Show, Enum)

-- | Result of a deriviation.
data Result a = Result {reason    :: !Reason
                       ,iters     :: !Int
                       ,lastDelta :: !Delta
                       ,lastCost  :: !Double
                       ,result    :: !a}
                deriving (Eq, Read, Show)

instance NFData a => NFData (Result a) where
    rnf = rnf . result




-- | A Dirichlet density.
data DirichletDensity = DD !(U.Vector Double)
                        deriving (Eq, Show)

instance NFData DirichletDensity where
    rnf = rwhnf

-- | @emptyDD n x@ is an \"empty\" Dirichlet density with size
--   @n@ and all alphas set to @x@.
emptyDD :: Int -> Double -> DirichletDensity
emptyDD = (DD .) . U.replicate

-- | Derive a Dirichlet density using a maximum likelihood method
--   as described by Karplus et al.  All training vectors should
--   have the same length, however this is not verified.
deriveDD :: DirichletDensity -> Predicate -> StepSize
         -> V.Vector TrainingVector -> Result DirichletDensity
deriveDD _ _ _ t | V.length t == 0 = error "Dirichlet.deriveDD: empty training data"
deriveDD (DD initial) (Pred maxIter' minDelta') (Step step) trainingData = train
    where
      size = U.length initial
      !trainingSize  = fromIntegral $ V.length trainingData
      trainingCounts = V.map (\t -> (t, U.sum t)) trainingData

      train =
        -- Initialization
        let as    = initial
            ws    = U.map log as
            sumAs = U.sum as
        in train' 0 1e100 sumAs (psi sumAs) ws as
      train' !iter !cost !sumAs !psiSumAs !ws !as =
        -- Reestimate alpha's
        let ws'    = U.izipWith calculateAlphas ws as
            as'    = U.map exp ws'
            sumAs' = U.sum as'
            cost'  = costDD' size as' sumAs' trainingCounts
            delta  = abs (cost' - cost)

        -- Verify convergence
        in case (delta <= minDelta', iter >= maxIter') of
             (True, _) -> Result Delta   iter delta cost' (DD as')
             (_, True) -> Result MaxIter iter delta cost' (DD as')
             _         -> train' (iter+1) cost' sumAs' (psi sumAs') ws' as'
       where
         calculateAlphas i w_old a_old =
           let s1 = trainingSize * (psiSumAs - psi a_old)
               f (t, sumT) = psi (t U.! i + a_old) - psi (sumT + sumAs)
           in w_old + step * a_old * (s1 + V.sum (V.map f trainingCounts))

-- | Cost function for deriving a Dirichlet density.  This
--   function is minimized by 'deriveDD'.
costDD :: DirichletDensity -> V.Vector TrainingVector -> Double
costDD (DD arr) tv = costDD' (U.length arr) arr (U.sum arr) $
                     V.map (\t -> (t, U.sum t)) tv

-- | 'costDD' needs to calculate the sum of all training vectors.
--   This functios avoids recalculting this quantity in
--   'deriveDD' multiple times.  This is the used by both
--   'costDD' and 'deriveDD'.
costDD' :: Int -> U.Vector Double -> Double -> V.Vector (TrainingVector, Double) -> Double
costDD' size alphas sumAs trainingCounts =
    let lngammaSumAs = lngamma sumAs
        f (t, sumT) = s + (U.sum $ U.zipWith g t alphas)
            where s = lngamma (sumT+1) + lngammaSumAs - lngamma (sumT + sumAs)
                  g t_i a_i = lngamma (t_i + a_i) - lngamma (t_i + 1) - lngamma a_i
    in negate . V.sum $ V.map f trainingCounts
{-# INLINE costDD' #-}