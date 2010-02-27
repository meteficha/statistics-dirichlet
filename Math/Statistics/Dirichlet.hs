---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet
-- Copyright   : (c) 2009 Felipe Lessa
-- License     : GPL
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
-- This module automatically uses
-- "Math.Statistics.Dirichlet.Density" or
-- "Math.Statistics.Dirichlet.Mixture" with the same API.
--
--------------------------------------------------------------------------

module Math.Statistics.Dirichlet
    (-- * Data types (re-exported)
     DirichletMixture(..)
    ,M.empty
    ,M.Component
    ,M.fromList
    ,M.toList
     -- * Options (re-exported)
    ,TrainingVector
    ,TrainingVectors
    ,StepSize(..)
    ,Delta
    ,Predicate(..)
    ,Reason(..)
    ,Result(..)
    -- * Functions
    ,derive
    ,cost
    -- * Util
    ,DirichletDensity(..)
    ,apply)
    where

import qualified Math.Statistics.Dirichlet.Density as D
import qualified Math.Statistics.Dirichlet.Mixture as M
import Math.Statistics.Dirichlet.Density (DirichletDensity(..))
import Math.Statistics.Dirichlet.Mixture (DirichletMixture(..))
import Math.Statistics.Dirichlet.Options



-- | Applies a function to a Dirichlet mixture.  If the mixture
-- has only one component, then a special function for densities
-- is applied.
apply :: (DirichletDensity -> a) -> (DirichletMixture -> a)
      -> DirichletMixture -> a
apply fd fm dm@(DM _ _ as)
      | M.dmComponents dm == 1 = fd (DD as)
      | otherwise              = fm dm


-- | Cost function for deriving a Dirichlet mixture (equation
-- 18).  This function is minimized by 'derive'.
cost :: TrainingVectors -> DirichletMixture -> Double
cost ns = apply (D.cost ns) (M.cost ns)

-- | Derive a Dirichlet mixture using a maximum likelihood method
-- as described by Karplus et al (equation 25).  All training
-- vectors should have the same length, however this is not
-- verified.
derive :: DirichletMixture -- ^ Initial mixture.
       -> Predicate        -- ^ Convergence test.
       -> StepSize -> TrainingVectors -> Result DirichletMixture
derive dm p s t = apply (\x -> change $ D.derive x p s t)
                        (\x ->          M.derive x p s t) dm
    where
      change :: Result DirichletDensity -> Result DirichletMixture
      change r = r {result = M.fromDD $ result r}
