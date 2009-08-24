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

module Math.Statistics.Dirichlet where

import Control.Applicative
import Control.Monad
import Control.Monad.ST
import Data.Array.Vector
import Numeric.GSL.Special.Gamma (lngamma, lnbeta)
import Numeric.GSL.Special.Psi (psi)



-- | Internal data type used to pair two doubles.  We don't use
--   the polymorphic (:*:) because GHC 6.10 can't unbox function
--   returns (see Don's statistics-fusion).
data DoublePair = (:**:) {-# UNPACK #-} !Double {-# UNPACK #-} !Double

-- | Logarithm of the beta function applied to a vector.
logBeta :: UArr Double -> Double
logBeta xs | lengthU xs == 2 = lnbeta (headU xs) (lastU xs)
           | otherwise       = fin $ foldlU go (0 :**: 0) xs
    where
      go  (val :**: sumXs) x = (val + lngamma x) :**: (sumXs + x)
      fin (val :**: sumXS)   = val - lngamma sumXS

-- | A vector used for deriving the parameters of a Dirichlet
--   density or mixture.
type TrainingVector = UArr Double

-- | Usually denoted by lowercase greek letter eta (η), size of
--   each step in the gradient. Should be greater than zero and
--   much less than one.
newtype StepSize = Step Double

-- | Predicate specifying when the training should be over.
data Predicate = Pred
    {maxIter  :: !Int    -- ^ Maximum number of iterations.
    ,minDelta :: !Double -- ^ Minimum delta to continue iterating.
    }

-- | Reason why the derivation was over.
data Reason = Delta | MaxIter

-- | Result of a deriviation.
data Result a = Result {reason    :: !Reason
                       ,iters     :: !Int
                       ,lastDelta :: !Double
                       ,result    :: !a}




-- | A Dirichlet density.
newtype DirichletDensity = DD (UArr Double)

-- | Returns a new filled 'MUArr'.
filled :: Int -> Double -> ST s (MUArr Double s)
filled size !value = newMU size >>= fill
    where fill !arr = go 0
              where go n | n == size = return arr
                         | otherwise = writeMU arr n value >> go (n+1)

-- | Sum of a mutable array.
sumMU :: Int -> MUArr Double s -> ST s Double
sumMU size !arr = go 0 0
    where go !x !acc | x == size = return acc
                     | otherwise = readMU arr x >>= go (x+1) . (acc+)

-- | Derive a Dirichlet density using a maximum likelihood method
--   as described by Karplus et al.  All training vectors should
--   have the same length, however this is not verified.
deriveDD :: Predicate -> StepSize
         -> [TrainingVector] -> ST s (Result DirichletDensity)
deriveDD _ _ [] = error "Dirichlet.deriveDD: empty training data"
deriveDD (Pred maxIter' minDelta') (Step step) trainingData
 = join $ return train `ap` arr 1 `ap` arr (exp 1) `ap` arr (exp 1)
    where
      arr  = filled size
      size = lengthU (head trainingData)
      !trainingSize  = fromIntegral $ length trainingData
      trainingCounts = map (\t -> t :*: sumU t) trainingData

      train = let sumAs = fromIntegral size * exp 1
              in train' 0 sumAs (psi sumAs)
      train' !iter !sumAs !psiSumAs !ws !as !new_as = do
              -- Reestimate alpha's
              (newDelta, newSumAs) <- calculateAlphas 0 0 0

              -- Verify convergence
              case (newDelta < minDelta', iter >= maxIter') of
                (True, _) -> Result Delta   iter newDelta <$> finish
                (_, True) -> Result MaxIter iter newDelta <$> finish
                _         -> train' (iter+1) newSumAs (psi newSumAs) ws new_as as
                             -- yes! we swap as with new_as
          where
            finish = DD <$> unsafeFreezeAllMU as

            calculateAlphas !i !newDelta !newSumAs
                | i == size = return (newDelta, newSumAs)
                | otherwise = do
              -- Old values
              w_old <- readMU ws i
              a_old <- readMU as i

              -- New values
              let s1 = trainingSize * (psiSumAs - psi a_old)
                  s2 = sum [psi (indexU t i + a_old) - psi (sumT + sumAs)
                                | t :*: sumT <- trainingCounts]
                  w_new = w_old + step * a_old * (s1 + s2)
                  a_new = exp w_new

              -- Write, calculate delta & sum
              writeMU ws     i w_new
              writeMU new_as i a_new
              let delta = abs (a_old - a_new)
              calculateAlphas (i+1) (max newDelta delta) (newSumAs + a_new)

