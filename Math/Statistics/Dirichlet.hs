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
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as MG
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as MU

import Control.Applicative
import Control.Monad
import Control.Monad.ST
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


--------------
sumMU :: Int -> MU.MVector s Double -> ST s Double
sumMU size vector = go size 0
    where
      go i s | s `seq` i < 0 = return s
             | otherwise     = do x <- MU.read vector i
                                  go i (s+x)

mutable :: U.Vector Double -> ST s (MU.MVector s Double)
mutable = MG.unstream . G.stream

immutable :: MU.MVector s Double -> ST s (U.Vector Double)
immutable v = do let n = MU.length v
                 r <- MU.new n
                 MU.copy r v
                 G.unsafeFreeze r
--------------



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
deriveDD (DD initial) (Pred maxIter' minDelta') (Step step) trainingData = runST train
    where
      size = U.length initial
      !trainingSize  = fromIntegral $ V.length trainingData
      trainingCounts = V.map (\t -> (t, U.sum t)) trainingData

      train = do
        -- Initialization
        as <- mutable initial
        ws <- mutable $ U.map log initial
        new_as <- MU.new size
        sumAs <- sumMU size as
        train' 0 1e100 sumAs (psi sumAs) ws as new_as
      train' :: Int -> Double -> Double -> Double -> MU.MVector s Double
             -> MU.MVector s Double -> MU.MVector s Double
             -> ST s (Result DirichletDensity)
      train' !iter !oldCost !sumAs !psiSumAs !ws !as !new_as = do
        -- Reestimate alpha's
        mapM_ calculateAlphas [0..size-1]
        newSumAs <- sumMU size new_as
        new_as'  <- immutable new_as
        let newCost = costDD' size new_as' newSumAs trainingCounts

        -- Verify convergence
        let delta = abs (newCost - oldCost)
        case (delta <= minDelta', iter >= maxIter') of
          (True, _) -> Result Delta   iter delta newCost <$> finish
          (_, True) -> Result MaxIter iter delta newCost <$> finish
          _         -> train' (iter+1) newCost newSumAs (psi newSumAs) ws new_as as
                       -- yes! we swap as with new_as
       where
         finish = DD <$> G.unsafeFreeze new_as
         calculateAlphas !i = do
           -- Old values
           w_old <- MU.read ws i
           a_old <- MU.read as i

           -- New values
           let s1 = trainingSize * (psiSumAs - psi a_old)
               s2 = V.map (\(t, sumT) -> psi (t U.! i + a_old) -
                                         psi (sumT + sumAs))
                          trainingCounts
               w_new = w_old + step * a_old * (s1 + V.sum s2)
           MU.write ws     i w_new
           MU.write new_as i (exp w_new) -- do not write in 'as'!

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