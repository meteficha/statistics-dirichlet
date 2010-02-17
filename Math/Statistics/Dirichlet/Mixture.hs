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
    (-- * Data types
     DirichletMixture(..)
    ,empty
    ,Component
    ,fromList
    ,toList
    -- * Functions used
    ,prob_a_n_theta)
    where

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Control.Parallel.Strategies (NFData(..), rwhnf)
import Numeric.GSL.Special.Gamma (lngamma, lnbeta)
import Numeric.GSL.Special.Psi (psi)

import qualified Math.Statistics.Dirichlet.Density as D
import Math.Statistics.Dirichlet.Density (DirichletDensity(..))
import Math.Statistics.Dirichlet.Options
import Math.Statistics.Dirichlet.Util



-- | A Dirichlet mixture.
data DirichletMixture =
    DM {dmWeights   :: {-# UNPACK #-} !(U.Vector Double)
       ,dmDensities :: {-# UNPACK #-} !(V.Vector DirichletDensity)}
        deriving (Eq)


instance Show DirichletMixture where
    showsPrec prec dm =
      showParen (prec > 10) $
      showString "fromList " .
      showsPrec 11 (toList dm)

instance Read DirichletMixture where
    readsPrec p ('(':xs) = let (ys,')':zs) = break (== ')') xs
                           in map (\(x,s) -> (x,s++zs)) $
                              readsPrec p ys
    readsPrec p xs = let [("fromList",list)] = lex xs
                     in map (\(x,s) -> (fromList x,s)) $
                        readsPrec p list

instance NFData DirichletMixture where
    rnf = rwhnf

-- | @empty q n x@ is an \"empty\" Dirichlet mixture with @q@
-- components.  Each component has size @n@, weight @1/q@ and all
-- alphas set to @x@.
empty :: Int -> Int -> Double -> DirichletMixture
empty q n x = let dd = D.empty n x
                  qs = recip $ fromIntegral q
              in DM {dmWeights   = U.replicate q qs
                    ,dmDensities = V.replicate q dd}
{-# INLINE empty #-}


-- | A list representation of a component of a Dirichlet mixture.
-- Used by 'fromList' and 'toList' only.
type Component = (Double, [Double])

-- | @fromList xs@ constructs a Dirichlet mixture from a
-- non-empty list of components.  Each component has a weight and
-- a list of alpha values.  The weights sum to 1, all lists must
-- have the same number of values and every number must be
-- non-negative.  All of these preconditions are verified for
-- clear mistakes.
fromList :: [Component] -> DirichletMixture
fromList components =
  let -- Vectors
      qs = U.fromList $ map               fst  components
      ds = V.fromList $ map (D.fromList . snd) components

      -- Properties of the mixture
      q  = length components
      n  = length (snd $ head components)

      -- Checks
      c0 = q >= 1
      c1 = abs (U.sum qs - 1) < 1e-2 -- we're quite permissive here
      c2 = U.all (>= 0) qs
      c3 = all ((== n) . length . snd) components
      c4 = all (all (>= 0)      . snd) components
      e  = error . ("Dirichlet.Mixture.fromList: " ++)
  in case (c0, c1, c2, c3, c4) of
       (True,_,_,_,_) -> e "there must be at least one component"
       (_,True,_,_,_) -> e "the sum of the weights must be one"
       (_,_,True,_,_) -> e "all weights must be greater than or equal to zero"
       (_,_,_,True,_) -> e "every component must have the same size"
       (_,_,_,_,True) -> e "all alphas must be greater than or equal to zero"
       _              -> DM qs ds

-- | @toList dm@ is the inverse of @fromList@, constructs a list
-- of components from a Dirichlet mixture.  There are no error
-- conditions and @toList . fromList == id@.
toList :: DirichletMixture -> [Component]
toList (DM qs ds) =
    let qs' = U.toList qs
        ds' = V.toList $ V.map D.toList ds
    in zip qs' ds'




-- | /Prob(a_j | n, theta)/ Defined in equation (16), "the
-- posterior probability of the /j/-th component of the mixture
-- given the vector of counts /n/".  We return the probabilities
-- for all /j/ in each vector.
--
-- Calculated as per equation (39) using 'logBeta'.
prob_a_n_theta :: TrainingVectors -> DirichletMixture -> V.Vector (U.Vector Double)
prob_a_n_theta ns (DM qs ds) =
    let -- Precalculate logBeta of all components
        !logBetaAlphas  = G.unstream $ G.stream $ V.map (logBeta . unDD) ds

        -- Calculates the factors for one of the training vectors.
        calc n i q lb_a = let a = unDD (ds V.! i)
                          in q * exp (logBeta (U.zipWith (+) n a) - lb_a)
        factors n       = let fs = U.izipWith (calc n) qs logBetaAlphas
                              total = U.sum fs
                          in U.map (/ total) fs
    in V.map factors ns
