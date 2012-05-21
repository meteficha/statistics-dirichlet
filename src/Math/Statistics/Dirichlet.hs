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
    ,empty
    ,Component
    ,fromList
    ,toList
     -- * Options (re-exported)
    ,TrainingVector
    ,TrainingVectors
    ,StepSize(..)
    ,Delta
    ,Predicate(..)
    ,Reason(..)
    ,Result(..)
    -- * Training data (re-exported)
    ,TrainingData
    ,prepareTraining
    -- * Functions (re-exported)
    ,derive
    ,cost)
    where

import Math.Statistics.Dirichlet.Mixture
import Math.Statistics.Dirichlet.Options
