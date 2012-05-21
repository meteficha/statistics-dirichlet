---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet
-- Copyright   : (c) 2009-2012 Felipe Lessa
-- License     : GPL
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
-- This module re-exports functions from
-- "Math.Statistics.Dirichlet.Mixture" and
-- "Math.Statistics.Dirichlet.Options" in a more digestable way.
-- Since this library is under-documented, I recommend reading
-- the documentation of the symbols re-exported here.
--
-- This module does not use "Math.Statistics.Dirichlet.Density"
-- in any way.  If you don't need mixtures then you should
-- probably use that module directly since it's faster and more
-- reliable (less magic happens there).
--
--------------------------------------------------------------------------

module Math.Statistics.Dirichlet
    ( -- * Data types (re-exported)
      DirichletMixture(..)
    , empty
    , Component
    , fromList
    , toList
      -- * Options (re-exported)
    , TrainingVector
    , TrainingVectors
    , StepSize(..)
    , Delta
    , Predicate(..)
    , Reason(..)
    , Result(..)
      -- * Training data (re-exported)
    , TrainingData
    , prepareTraining
      -- * Functions (re-exported)
    , derive
    , cost
    ) where

import Math.Statistics.Dirichlet.Mixture
import Math.Statistics.Dirichlet.Options
