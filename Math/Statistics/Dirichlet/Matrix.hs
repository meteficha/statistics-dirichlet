---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet.Matrix
-- Copyright   : (c) 2009 Felipe Lessa
-- License     : GPL
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
-- Implement matrices using plain 'U.Vector's with data stored in
-- row-major order (i.e. the first elements correspond to the
-- first row).
--
--------------------------------------------------------------------------

module Math.Statistics.Dirichlet.Matrix
    (-- * Data type
     Matrix(..)
    ,size
     -- * Constructing
    ,replicate
    ,replicateRows
    ,fromVector
     -- * Rows
    ,rows
    ,(!!!)
    ,rowMap
     -- * Columns
    ,cols
    ,col
     -- * Maps and zips
    ,umap
    ,map
    ,imap
    ,uzipWith
    ,zipWith
    ,izipWith
    ) where

import Prelude hiding (replicate, map, zipWith)
import qualified Data.Vector.Fusion.Stream as S
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U


-- | A matrix.
data Matrix = M {mRows :: {-# UNPACK #-} !Int
                ,mCols :: {-# UNPACK #-} !Int
                ,mData :: {-# UNPACK #-} !(U.Vector Double)}
            deriving (Eq, Ord)

-- | Size of the matrix.
size :: Matrix -> (Int,Int)
size m = (mRows m, mCols m)



-- | A matrix where all elements are of the same value.
replicate :: (Int,Int) -> Double -> Matrix
replicate (r,c) v = M {mRows = r
                      ,mCols = c
                      ,mData = U.replicate (r*c) v}

-- | A matrix where all rows are of the same value.
replicateRows :: Int -> U.Vector Double -> Matrix
replicateRows r v =
    let c = U.length v
    in M {mRows = r
         ,mCols = c
         ,mData = U.generate (r*c) (\i -> v U.! (i `mod` c))}

-- | Creates a matrix from a vector of vectors.  It is not
-- verified that the vectors have the right length.
fromVector :: (G.Vector v (w Double), G.Vector w Double)
           => (Int,Int) -> v (w Double) -> Matrix
fromVector (r,c) v =
    M {mRows = r
      ,mCols = c
      ,mData = G.unstream $ S.concatMap G.stream $ G.stream v}
{-# INLINE fromVector #-}




-- | /O(rows)/ Rows of the matrix.  Each element takes /O(1)/ time and
-- storage.
rows :: Matrix -> [U.Vector Double]
rows m = do
  i <- [0, mCols m .. U.length (mData m) - mCols m]
  return $ U.unsafeSlice i (mCols m) (mData m)

-- | /O(1)/ @m !!! i@ is the @i@-th row of the matrix.
(!!!) :: Matrix -> Int -> U.Vector Double
m !!! i = U.slice (i * mCols m) (mCols m) (mData m)

rowMap :: (U.Vector Double -> Double) -> Matrix -> U.Vector Double
rowMap f m = U.generate (mRows m) (f . s)
    where s i = U.unsafeSlice (i * mCols m) (mCols m) (mData m)





-- | /O(rows*cols)/ Columns of the matrix.  Each element takes
-- /O(rows)/ time and storage.
cols :: Matrix -> [U.Vector Double]
cols m = [m `col` i | i <- [0..mCols m - 1]]

-- | /O(rows)/ @m `col` i@ is the @i@-th column of the matrix.
col :: Matrix -> Int -> U.Vector Double
m `col` i = U.backpermute (mData m) $ U.enumFromStepN i (mCols m) (mRows m)






umap :: (U.Vector Double -> U.Vector Double) -> Matrix -> Matrix
umap f m = m {mData = f (mData m)}

map :: (Double -> Double) -> Matrix -> Matrix
map f = umap (U.map f)

imap :: ((Int,Int) -> Double -> Double) -> Matrix -> Matrix
imap f m = umap (U.imap (f . indices m)) m

uzipWith :: (U.Vector Double -> U.Vector Double -> U.Vector Double)
         -> Matrix -> Matrix -> Matrix
uzipWith f m n
    | mRows m /= mRows n = materror "uzipWith" "mRows"
    | mCols m /= mCols n = materror "uzipWith" "mCols"
    | otherwise          = m {mData = f (mData m) (mData n)}

zipWith :: (Double -> Double -> Double) -> Matrix -> Matrix -> Matrix
zipWith f = uzipWith (U.zipWith f)

izipWith :: ((Int,Int) -> Double -> Double -> Double)
         -> Matrix -> Matrix -> Matrix
izipWith f m = uzipWith (U.izipWith (f . indices m)) m

indices :: Matrix -> Int -> (Int, Int)
indices m i = i `divMod` mCols m



materror :: String -> String -> a
materror f e = error $ "Math.Statistics.Dirichlet.Matrix." ++ f ++ ": " ++ e