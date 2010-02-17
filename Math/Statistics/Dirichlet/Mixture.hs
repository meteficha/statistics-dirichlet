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
    (DirichletDensity(..)
    ,emptyDD
    ,deriveDD
    ,costDD

    ,logBeta)
    where

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Control.Parallel.Strategies (NFData(..), rwhnf)
import Numeric.GSL.Special.Gamma (lngamma, lnbeta)
import Numeric.GSL.Special.Psi (psi)

import Math.Statistics.Dirichlet.Options



-- | Logarithm of the beta function applied to a vector.
logBeta :: U.Vector Double -> Double
logBeta xs | U.length xs == 2 = lnbeta (U.head xs) (U.last xs)
           | otherwise        = U.sum (U.map lngamma xs) - lngamma (U.sum xs)
{-# INLINE logBeta #-}

-- | A Dirichlet density.
newtype DirichletDensity = DD (U.Vector Double) deriving (Eq)

instance Show DirichletDensity where
    showsPrec prec (DD v) =
      showParen (prec > 10) $
      showString "listDD " .
      showsPrec 11 (U.toList v)

instance Read DirichletDensity where
    readsPrec p ('(':xs) = let (ys,')':zs) = break (== ')') xs
                           in map (\(x,s) -> (x,s++zs)) $
                              readsPrec p ys
    readsPrec p xs = let [("listDD",list)] = lex xs
                     in map (\(x,s) -> (listDD x,s)) $
                        readsPrec p list

instance NFData DirichletDensity where
    rnf = rwhnf

-- | @emptyDD n x@ is an \"empty\" Dirichlet density with size
--   @n@ and all alphas set to @x@.
emptyDD :: Int -> Double -> DirichletDensity
emptyDD = (DD .) . U.replicate
{-# INLINE emptyDD #-}

-- | @listDD xs@ constructs a Dirichlet density from a list of
-- alpha values.
listDD :: [Double] -> DirichletDensity
listDD = DD . U.fromList

infinity :: Double
infinity = 1e100

-- | Derive a Dirichlet density using a maximum likelihood method
--   as described by Karplus et al.  All training vectors should
--   have the same length, however this is not verified.
deriveDD :: DirichletDensity -> Predicate -> StepSize
         -> TrainingVectors -> Result DirichletDensity
deriveDD (DD initial) (Pred maxIter' minDelta_ deltaSteps')
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
      err = error . ("Dirichlet.deriveDD: " ++)

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
      train = train' 0 infinity (U.sum initial) $
              U.map (\x -> (log x, x)) initial

      train' !iter !cost !sumAs !alphas =
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
            !cost'    = if calcCost then costDD' (snd $ U.unzip alphas') sumAs'
                                                 trainingData trainingSums
                                    else cost -- use old cost
            !delta    = abs (cost' - cost)

        -- Verify convergence.  Even with MaxIter we only stop
        -- iterating if the delta was calculated.  Otherwise we
        -- wouldn't be able to tell the caller why the delta was
        -- still big when we reached the limit.
        in case (calcCost, delta <= minDelta', iter >= maxIter') of
             (True, True, _) -> Result Delta   iter delta cost' $ finish alphas'
             (True, _, True) -> Result MaxIter iter delta cost' $ finish alphas'
             _               -> train' (iter+1) cost' sumAs' alphas'

-- | Cost function for deriving a Dirichlet density.  This
--   function is minimized by 'deriveDD'.
costDD :: DirichletDensity -> TrainingVectors -> Double
costDD (DD arr) tv = costDD' arr (U.sum arr) tv $
                     G.unstream $ G.stream $ V.map U.sum tv

-- | 'costDD' needs to calculate the sum of all training vectors.
--   This functios avoids recalculting this quantity in
--   'deriveDD' multiple times.  This is the used by both
--   'costDD' and 'deriveDD'.
costDD' :: U.Vector Double -> Double -> TrainingVectors -> U.Vector Double -> Double
costDD' !alphas !sumAs !trainingData !trainingSums =
    let !lngammaSumAs = lngamma sumAs
        f t = U.sum $ U.zipWith w t alphas
            where w t_i a_i = lngamma (t_i + a_i) - lngamma (t_i + 1) - lngamma a_i
        g sumT = lngamma (sumT+1) + lngammaSumAs - lngamma (sumT + sumAs)
    in negate $ (V.sum $ V.map f trainingData)
              + (U.sum $ U.map g trainingSums)
{-# INLINE costDD' #-}