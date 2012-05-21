---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet.Util
-- Copyright   : (c) 2009 Felipe Lessa
-- License     : GPL
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
--------------------------------------------------------------------------

module Math.Statistics.Dirichlet.Util
    (infinity
    ,logBeta)
    where

import qualified Data.Vector.Unboxed as U
import Numeric.GSL.Special.Gamma (lngamma, lnbeta)



-- | Logarithm of the beta function applied to a vector.
logBeta :: U.Vector Double -> Double
logBeta xs | U.length xs == 2 = lnbeta (U.head xs) (U.last xs)
           | otherwise        = U.sum (U.map lngamma xs) - lngamma (U.sum xs)

-- | Infinity, currently defined as @1e100@.  Used mainly as the
-- initial cost.
infinity :: Double
infinity = 1e100
{-# INLINE infinity #-}