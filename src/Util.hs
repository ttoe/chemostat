{-# LANGUAGE NoImplicitPrelude #-}

module Util
  ( diff
  , getMatrixCol
  , zipVecWith2
  , matrixAddSumCol
  , times
  ) where

import Foundation
import Foundation.Numerical (roundDown)

import Data.List ((!!), (++))

import qualified Data.Vector.Storable as VS

import qualified Numeric.LinearAlgebra as LA
import Numeric.GSL.Differentiation (derivCentral)

type D = Double

-- alias differentiation function to central derivative with initial step size 0.01
diff :: (D -> D) -> D -> D
diff fun point = fst $ derivCentral 0.01 fun point

-- takes the time vector and the solution matrix
-- creates a list of lists containing tuples
-- for each variable and the corresponding time
-- necessary to satisfy plotting funtion
-- makePlottableTuples :: DV.Vector Double -> LA.Matrix Double -> VS.Vector (Double, Double)
-- makePlottableTuples vec mtx = fmap (zip (LA.toList vec) . LA.toList) (toColumns mtx)
-- makePlottableTuples vec mtx = DV.map (zipWith (,) vec) (LA.toColumns mtx)

-- create a string like "2017-06-09_131211"
-- as a timestamp to use for writing files
-- getNowTimeString :: IO ()
-- getNowTimeString = do
--   now <- getCurrentTime
--   pure (formatTime defaultTimeLocale "%F_%H%M%S" now)
--
-- -- zipColumns :: (V.Storable a, V.Storable b, LA.Element a) =>
--   -- (a -> a -> b) -> LA.Matrix a -> Int -> Int -> V.Vector b
-- zipColumns :: (Double -> Double -> Double) -> LA.Matrix Double -> Int -> Int -> VS.Vector Double
-- zipColumns f m x y = DV.zipWith f (columns !! x) (columns !! y)
--   where
--     columns = LA.toColumns m
--
-- -- zipWithColumn :: (V.Storable a, V.Storable b, V.Storable c, LA.Element a) =>
--   -- (a -> b -> c) -> LA.Matrix a -> V.Vector b -> Int -> V.Vector c
-- zipWithColumn :: (Double -> Double -> Double) -> LA.Matrix Double -> LA.Vector Double -> Int -> LA.Vector Double
-- zipWithColumn fun mtx vec col = DV.zip vec mcol
--   where
--     mcol = (LA.toColumns mtx) !! col

zipVecWith2 :: (D -> D -> D) -> LA.Vector D -> LA.Vector D -> LA.Vector D
zipVecWith2 f x y = VS.zipWith f x y

getMatrixCol :: Int -> LA.Matrix D -> LA.Vector D
getMatrixCol i m = (LA.toColumns m) !! i

matrixAddSumCol :: Int -> Int -> LA.Matrix D -> LA.Matrix D
matrixAddSumCol x y m = LA.fromColumns (columns ++ [colSum])
  where
    columns = LA.toColumns m
    colSum  = VS.zipWith (+) (columns !! x) (columns !! y)

times :: D -> D -> D -> LA.Vector D
times from to stepsize = LA.linspace steps (from, to)
  where
    steps = roundDown $ (to-from)/stepsize
