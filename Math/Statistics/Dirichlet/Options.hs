---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet.Options
-- Copyright   : (c) 2009 Felipe Lessa
-- License     : GPL
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
--------------------------------------------------------------------------

module Math.Statistics.Dirichlet.Options
    (TrainingVector
    ,TrainingVectors
    ,StepSize(..)
    ,Delta
    ,Predicate(..)
    ,Reason(..)
    ,Result(..))
    where

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

-- | A vector used for deriving the parameters of a Dirichlet
--   density or mixture.
type TrainingVector = U.Vector Double

-- | A vector of training vectors.  This is the only vector that
-- is not unboxed (for obvious reasons).
type TrainingVectors = V.Vector TrainingVector

-- | Usually denoted by lowercase greek letter eta (η), size of
--   each step in the gradient. Should be greater than zero and
--   much less than one.
newtype StepSize = Step Double

-- | Maximum difference between costs to consider that the
--   process converged.
type Delta = Double

-- | Predicate specifying when the training should be over.
data Predicate = Pred
    {maxIter    :: !Int    -- ^ Maximum number of iterations.
    ,minDelta   :: !Delta  -- ^ Minimum delta to continue iterating.
                           --   This is invariant of @deltaSteps@, which
                           --   means that if @deltaSteps@ is @2@ then
                           --   minDelta will be considered twice bigger
                           --   to account for the different @deltaSteps@.
    ,deltaSteps :: !Int    -- ^ How many estimation steps should be done
                           --   before recalculating the delta.  If
                           --   @deltaSteps@ is @1@ then it will be
                           --   recalculated on every step.
    }
                 deriving (Eq, Read, Show)

-- | Reason why the derivation was over.
data Reason = Delta    -- ^ The difference between applications
                       -- of the cost function dropped below the
                       -- minimum delta.  In other words, it
                       -- coverged.
            | MaxIter  -- ^ The maximum number of iterations was
                       -- reached while the delta was still
                       -- greater than the minimum delta.
              deriving (Eq, Read, Show, Enum)

-- | Result of a deriviation.
data Result a =
    Result {reason    :: !Reason  -- ^ Reason why the derivation was over.
           ,iters     :: !Int     -- ^ Number of iterations spent.
           ,lastDelta :: !Delta   -- ^ Last difference between costs.
           ,lastCost  :: !Double  -- ^ Last cost (i.e. the cost of the result).
           ,result    :: !a       -- ^ Result obtained.
           }
    deriving (Eq, Read, Show)

instance NFData a => NFData (Result a) where
    rnf = rnf . result
