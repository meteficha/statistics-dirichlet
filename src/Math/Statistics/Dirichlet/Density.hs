---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet.Density
-- Copyright   : (c) 2009-2012 Felipe Lessa
-- License     : BSD3
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
--------------------------------------------------------------------------

module Math.Statistics.Dirichlet.Density
    ( DirichletDensity(..)
    , empty
    , fromList
    , toList
    , derive
    , cost
    ) where

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Control.DeepSeq (NFData(..))
import Numeric.GSL.Special.Gamma (lngamma)
import Numeric.GSL.Special.Psi (psi)

import Math.Statistics.Dirichlet.Options
import Math.Statistics.Dirichlet.Util



-- | A Dirichlet density.
newtype DirichletDensity = DD {unDD :: U.Vector Double} deriving (Eq)

instance Show DirichletDensity where
    showsPrec prec (DD v) =
      showParen (prec > 10) $
      showString "fromList " .
      showsPrec 11 (U.toList v)

instance Read DirichletDensity where
    readsPrec p ('(':xs) = let (ys,')':zs) = break (== ')') xs
                           in map (\(x,s) -> (x,s++zs)) $
                              readsPrec p ys
    readsPrec p xs = let [("fromList",list)] = lex xs
                     in map (\(x,s) -> (fromList x,s)) $
                        readsPrec p list

instance NFData DirichletDensity where
    rnf DD {} = ()

-- | @empty n x@ is an \"empty\" Dirichlet density with size
-- @n@ and all alphas set to @x@.
empty :: Int -> Double -> DirichletDensity
empty = (DD .) . U.replicate
{-# INLINE empty #-}

-- | @fromList xs@ constructs a Dirichlet density from a list of
-- alpha values.
fromList :: [Double] -> DirichletDensity
fromList = DD . U.fromList
{-# INLINE fromList #-}

-- | @toList d@ deconstructs a Dirichlet density to a list of
-- alpha values.
toList :: DirichletDensity -> [Double]
toList (DD xs) = U.toList xs
{-# INLINE toList #-}

-- | Derive a Dirichlet density using a maximum likelihood method
-- as described by Karplus et al (equation 26).  All training
-- vectors should have the same length, however this is not
-- verified.
derive :: DirichletDensity -> Predicate -> StepSize
         -> TrainingVectors -> Result DirichletDensity
derive (DD initial) (Pred maxIter' minDelta_ deltaSteps' _ _)
             (Step step) trainingData
    | V.length trainingData == 0 = err "empty training data"
    | U.length initial < 1       = err "empty initial vector"
    | maxIter' < 1               = err "non-positive maxIter"
    | minDelta_ < 0              = err "negative minDelta"
    | deltaSteps' < 1            = err "non-positive deltaSteps"
    | step <= 0                  = err "non-positive step"
    | step >= 1                  = err "step greater than one"
    | otherwise                  = train
    where
      err = error . ("Dirichlet.derive: " ++)

      -- Compensate the different deltaSteps.
      !minDelta'    = minDelta_ * fromIntegral deltaSteps'

      -- Number of training sequences.
      !trainingSize = fromIntegral $ V.length trainingData

      -- Sums of each training sequence.
      trainingSums :: U.Vector Double
      !trainingSums = G.unstream $ G.stream $ V.map U.sum trainingData

      -- Functions that work on the alphas only (and not their logs).
      calcSumAs = U.sum . snd . U.unzip
      finish    = DD    . snd . U.unzip

      -- Start training in the zero-th iteration and with
      -- infinite inital cost.
      train = train' 1 infinity (U.sum initial) $
              U.map (\x -> (log x, x)) initial

      train' !iter !oldCost !sumAs !alphas =
        -- Reestimate alpha's.
        let !alphas'  = U.imap calculateAlphas alphas
            !psiSumAs = psi sumAs
            !psiSums  = U.sum $ U.map (\sumT -> psi $ sumT + sumAs) trainingSums
            calculateAlphas !i (!w, !a) =
              let !s1 = trainingSize * (psiSumAs - psi a)
                  !s2 = V.sum $ V.map (\t -> psi $ t U.! i + a) trainingData
                  !w' = w + step * a * (s1 + s2 - psiSums)
                  !a' = exp w'
              in (w', a')

        -- Recalculate constants.
            !sumAs'   = calcSumAs alphas'
            !calcCost = iter `mod` deltaSteps' == 0
            !cost'    = if calcCost then newCost else oldCost
             where newCost = costWorker (snd $ U.unzip alphas') sumAs'
                                        trainingData trainingSums
            !delta    = abs (cost' - oldCost)

        -- Verify convergence.  Even with MaxIter we only stop
        -- iterating if the delta was calculated.  Otherwise we
        -- wouldn't be able to tell the caller why the delta was
        -- still big when we reached the limit.
        in case (calcCost, delta <= minDelta', iter >= maxIter') of
             (True, True, _) -> Result Delta   iter delta cost' $ finish alphas'
             (True, _, True) -> Result MaxIter iter delta cost' $ finish alphas'
             _               -> train' (iter+1) cost' sumAs' alphas'

-- | Cost function for deriving a Dirichlet density (equation
-- 18).  This function is minimized by 'derive'.
cost :: TrainingVectors -> DirichletDensity -> Double
cost tv (DD arr) = costWorker arr (U.sum arr) tv $
                     G.unstream $ G.stream $ V.map U.sum tv

-- | 'cost' needs to calculate the sum of all training vectors.
-- This functios avoids recalculting this quantity in 'derive'
-- multiple times.  This is the used by both 'cost' and 'derive'.
costWorker :: U.Vector Double -> Double -> TrainingVectors -> U.Vector Double -> Double
costWorker !alphas !sumAs !trainingData !trainingSums =
    let !lngammaSumAs = lngamma sumAs
        f t = U.sum $ U.zipWith w t alphas
            where w t_i a_i = lngamma (t_i + a_i) - lngamma (t_i + 1) - lngamma a_i
        g sumT = lngamma (sumT+1) - lngamma (sumT + sumAs)
    in negate $ (V.sum $ V.map f trainingData)
              + (U.sum $ U.map g trainingSums)
              + lngammaSumAs * fromIntegral (U.length trainingSums)
{-# INLINE costWorker #-}
