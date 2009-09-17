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
                 deriving (Eq, Show)

-- | Reason why the derivation was over.
data Reason = Delta | MaxIter
              deriving (Eq, Show, Enum)

-- | Result of a deriviation.
data Result a = Result {reason    :: !Reason
                       ,iters     :: !Int
                       ,lastDelta :: !Double
                       ,lastCost  :: !Double
                       ,result    :: !a}
                deriving (Eq, Show)





-- | Sum of a mutable array.  Uses Kahan summation algorithm.
kahanSumMU :: MUArr Double s -> ST s Double
kahanSumMU !arr | size >= 1 = readMU arr 0 >>= go 1 0
                | otherwise = return 0
    where go !i !comp !acc
              | i == size = return acc
              | otherwise = do x <- readMU arr i
                               let x' = x - comp
                                   acc'  = acc + x'
                                   comp' = (acc' - acc) - x'
                               go (i+1) comp' acc'
          size = lengthMU arr

-- | Sum of a list.  Uses Kahan summation algorithm.
kahanSum :: [Double] -> Double
kahanSum []     = 0
kahanSum (y:ys) = go 0 y ys
    where go !comp !acc (x:xs) = let x' = x - comp
                                     acc'  = acc + x'
                                     comp' = (acc' - acc) - x'
                                 in go comp' acc' xs
          go _ acc [] = acc





-- | A Dirichlet density.
data DirichletDensity = DD !(UArr Double)
                        deriving (Eq, Show)

-- | @emptyDD n x@ is an \"empty\" Dirichlet density with size
--   @n@ and all alphas set to @x@.
emptyDD :: Int -> Double -> DirichletDensity
emptyDD = (DD .) . replicateU

-- | Derive a Dirichlet density using a maximum likelihood method
--   as described by Karplus et al.  All training vectors should
--   have the same length, however this is not verified.
deriveDD :: DirichletDensity -> Predicate -> StepSize
         -> [TrainingVector] -> Result DirichletDensity
deriveDD _ _ _ [] = error "Dirichlet.deriveDD: empty training data"
deriveDD (DD initial) (Pred maxIter' minDelta') (Step step) trainingData = runST train
    where
      size = lengthU (head trainingData)
      !trainingSize  = fromIntegral $ length trainingData
      trainingCounts = map (\t -> t :*: sumU t) trainingData

      train = do
        -- Initialization
        [ws, as, new_as] <- replicateM 3 (newMU size)
        copyMU as 0 initial
        copyMU ws 0 (mapU log initial)
        sumAs <- kahanSumMU as
        train' 0 1e100 sumAs (psi sumAs) ws as new_as
      train' !iter !oldCost !sumAs !psiSumAs !ws !as !new_as = do
        -- Reestimate alpha's
        mapM_ calculateAlphas [0..size-1]
        newSumAs <- kahanSumMU new_as
        newCost  <- costDD' new_as newSumAs trainingCounts

        -- Verify convergence
        let delta = abs (newCost - oldCost)
        case (delta <= minDelta', iter >= maxIter') of
          (True, _) -> Result Delta   iter delta newCost <$> finish
          (_, True) -> Result MaxIter iter delta newCost <$> finish
          _         -> train' (iter+1) newCost newSumAs (psi newSumAs) ws new_as as
                       -- yes! we swap as with new_as
       where
         finish = DD <$> unsafeFreezeAllMU new_as
         calculateAlphas !i = do
           -- Old values
           w_old <- readMU ws i
           a_old <- readMU as i

           -- New values
           let s1 = trainingSize * (psiSumAs - psi a_old)
               s2 = do t :*: sumT <- trainingCounts
                       [psi (t `indexU` i + a_old), -psi (sumT + sumAs)]
               w_new = w_old + step * a_old * kahanSum (s1 : s2)
           writeMU ws     i w_new
           writeMU new_as i (exp w_new) -- do not write in 'as'!

-- | Cost function for deriving a Dirichlet density.  This function is minimized.
costDD :: DirichletDensity -> [TrainingVector] -> Double
costDD (DD arr) tv = runST (go $ map (\t -> t :*: sumU t) tv)
    where
      size = lengthU arr
      go tc = do alphas <- newMU size
                 copyMU alphas 0 arr
                 sumAs <- kahanSumMU alphas
                 costDD' alphas sumAs tc

costDD' :: MUArr Double s -> Double -> [TrainingVector :*: Double] -> ST s Double
costDD' alphas sumAs trainingCounts =
    let lngammaSumAs = lngamma sumAs
        size = lengthMU alphas
        f (t :*: sumT) = go 0 $ lngamma (sumT+1) + lngammaSumAs - lngamma (sumT + sumAs)
            where go !i !acc | i == size = return acc
                             | otherwise = do
                    a_i <- readMU alphas i
                    let t_i = t `indexU` i
                    go (i+1) $ acc + lngamma (t_i + a_i) - lngamma (t_i + 1) - lngamma a_i
    in (negate . kahanSum) `fmap` mapM f trainingCounts
