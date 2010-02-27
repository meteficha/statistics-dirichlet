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
    (-- * Basic
     Matrix(..)
    ,size
    ,(!)
     -- * Constructing
    ,replicate
    ,replicateRows
    ,fromVector
    ,fromVectorT
     -- * Rows
    ,rows
    ,(!!!)
     -- * Columns
    ,cols
    ,col
     -- * Maps and zips
    ,umap
    ,map
    ,imap
    ,rowmap
    ,uzipWith
    ,zipWith
    ,izipWith
    ,rzipWith
     -- * Other
    ,transpose
    ) where

import Prelude hiding (replicate, map, zipWith)
import System.IO.Unsafe (unsafePerformIO)
import qualified Data.Vector as V
import qualified Data.Vector.Fusion.Stream as S
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as MU


-- | A matrix.
data Matrix = M {mRows :: {-# UNPACK #-} !Int
                ,mCols :: {-# UNPACK #-} !Int
                ,mData :: {-# UNPACK #-} !(U.Vector Double)}
            deriving (Eq, Ord, Show)

-- | Size of the matrix.
size :: Matrix -> (Int,Int)
size m = (mRows m, mCols m)

-- | Element at position.
(!) :: Matrix -> (Int,Int) -> Double
(!) m (r,c) = mData m U.! (r * mCols m + c)



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

-- | Creates a matrix from a vector of vectors.  It *is not*
-- verified that the vectors have the right length.
fromVector :: (G.Vector v (w Double), G.Vector w Double)
           => v (w Double) -> Matrix
fromVector v =
    M {mRows = G.length v
      ,mCols = G.length (G.head v)
      ,mData = G.unstream $ S.concatMap G.stream $ G.stream v}

-- | Creates a matrix from a vector of vectors.  The vectors are
-- transposed, so @fromVectorT@ is the same as @transpose
-- . fromVector@. It *is* verified that the vectors have the
-- right length.
fromVectorT :: (G.Vector v (w Double), G.Vector w Double)
           => v (w Double) -> Matrix
fromVectorT v =
    M {mRows = c
      ,mCols = r
      ,mData = unsafePerformIO $ do
                 m <- MU.new (r*c)
                 fillCol m r
                 G.unsafeFreeze m}
  where
    r = G.length v
    c = G.length (G.head v)
    fillCol _ 0 = return ()
    fillCol m j = let j' = j-1
                  in fillRow m (v G.! j') j' c >> fillCol m j'
    fillRow _ _   _  0 = return ()
    fillRow m clm j' i = let i' = i-1
                             x  = clm G.! i'
                         in MU.write m (i' * r + j') x >> fillRow m clm j' i'




-- | /O(rows)/ Rows of the matrix.  Each element takes /O(1)/ time and
-- storage.
rows :: Matrix -> V.Vector (U.Vector Double)
rows m = G.map (\i -> U.unsafeSlice i (mCols m) (mData m)) $
         G.enumFromStepN 0 (mCols m) (mRows m)

-- | /O(1)/ @m !!! i@ is the @i@-th row of the matrix.
(!!!) :: Matrix -> Int -> U.Vector Double
m !!! i = U.slice (i * mCols m) (mCols m) (mData m)






-- | /O(rows*cols)/ Columns of the matrix.  Each element takes
-- /O(rows)/ time and storage.
cols :: Matrix -> V.Vector (U.Vector Double)
cols m = V.generate (mCols m) (m `col`)

-- | /O(rows)/ @m `col` i@ is the @i@-th column of the matrix.
col :: Matrix -> Int -> U.Vector Double
m `col` i = U.backpermute (mData m) $ U.enumFromStepN i (mCols m) (mRows m)






umap :: (U.Vector Double -> U.Vector Double) -> Matrix -> Matrix
umap f m = m {mData = f (mData m)}

map :: (Double -> Double) -> Matrix -> Matrix
map f = umap (U.map f)

imap :: ((Int,Int) -> Double -> Double) -> Matrix -> Matrix
imap f m = umap (U.imap (f . indices m)) m

rowmap :: (U.Vector Double -> Double) -> Matrix -> U.Vector Double
rowmap f m = U.generate (mRows m) (f . s)
    where s i = U.unsafeSlice (i * mCols m) (mCols m) (mData m)

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

-- | @rzipWith f m n@ is a matrix with the same number of rows as
-- @m@.  The @i@-th row is obtained by applying @f@ to the @i@-th
-- rows of @m@ and @n@.
rzipWith :: (Int -> U.Vector Double -> U.Vector Double -> U.Vector Double)
          -> Matrix -> Matrix -> Matrix
rzipWith f m n
    | mRows m /= mRows n = materror "rzipWithN" "mRows"
    | mCols m /= mCols n = materror "rzipWithN" "mCols"
    | otherwise          = fromVector $ V.izipWith f (rows m) (rows n)


indices :: Matrix -> Int -> (Int, Int)
indices m i = i `divMod` mCols m



transpose :: Matrix -> Matrix
transpose m =
    let f i = let (r,c) = i `divMod` mRows m
              in m ! (c,r)
    in M {mRows = mCols m
         ,mCols = mRows m
         ,mData = U.generate (mRows m * mCols m) f}

{-# RULES
  "transpose/transpose"   forall m. transpose (transpose m) = m;
  "transpose/fromVector"  forall v. transpose (fromVector v) = fromVectorT v;
  "transpose/fromVectorT" forall v. transpose (fromVectorT v) = fromVector v;
  #-}



materror :: String -> String -> a
materror f e = error $ "Math.Statistics.Dirichlet.Matrix." ++ f ++ ": " ++ e