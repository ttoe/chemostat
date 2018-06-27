module Util
  ( diff
  , matrixAddSumCol
  , times
  , nowTimeString
  , writePlot
  , mkPlotTuples
  , matrixAddTimeCol
  , mtxLast
  , printTimeDiff
  -- , findLocMaxWithIx
  -- , findLocMinWithIx
  , findLocMinMaxWithIx
  ) where

import Data.Time.Clock (UTCTime, getCurrentTime, diffUTCTime)
import Data.Time.Format (formatTime, defaultTimeLocale)
import qualified Data.Vector.Storable as VS
import Graphics.Rendering.Chart.Backend.Cairo (FileFormat(..), toFile, _fo_format)
import Graphics.Rendering.Chart.Easy (Default, ToRenderable, EC, def)
import Numeric.LinearAlgebra ( Vector
                             , Matrix
                             , Extractor(..)
                             , (??)
                             , linspace
                             , toColumns
                             , fromColumns
                             , toList
                             , atIndex
                             )
import Numeric.GSL.Differentiation (derivCentral)
import Data.List (zip4)

type D = Double


-- alias differentiation function to central derivative with initial step size 0.01
diff :: (D -> D) -> D -> D
diff fun point = fst $ derivCentral 0.01 fun point


-- takes a matrix where column 0 is the time
-- and the other columns are the solutions at a given time
-- creates a list of lists containing tuples
-- for each variable and the corresponding time
mkPlotTuples :: Matrix D -> [[(D, D)]]
mkPlotTuples m = fmap (zip timeCol) columns
  where
    columns = fmap toList $ toColumns m
    timeCol = columns !! 0


mtxLast :: Int -> Matrix D -> Matrix D
mtxLast l m = m ?? (TakeLast takeLast, All)
  where
    timeStep = m `atIndex` (1, 0) -- 2nd value (1, _) in 1st column (_, 0) is time step from times function
    takeLast = floor $ fromIntegral l / timeStep


-- create a timestamp like "2017-06-09_131211"
nowTimeString :: IO String
nowTimeString = do
  now <- getCurrentTime
  pure (formatTime defaultTimeLocale "%F_%H%M%S" now)


matrixAddSumCol :: Int -> Int -> Matrix D -> Matrix D
matrixAddSumCol x y m = fromColumns (columns <> [colSum])
  where
    columns = toColumns m
    colSum  = VS.zipWith (+) (columns !! x) (columns !! y)


matrixAddTimeCol :: Vector D -> Matrix D -> Matrix D
matrixAddTimeCol t m = fromColumns ([t] <> toColumns m)


times :: D -> D -> D -> Vector D
times from to stepsize = linspace steps (from, to)
  where
    steps = floor $ (to-from)/stepsize


printTimeDiff :: UTCTime -> UTCTime -> IO ()
printTimeDiff startTime endTime = print $ diffUTCTime endTime startTime


-- these 3 min/max finding functions could be unified with a selector

-- findLocMaxWithIx :: (Ord a, Num a) => [a] -> [(a, Int)]
-- findLocMaxWithIx xs = map keepValueAndIx . filter isLocalMaximum $ zippedToWindow
--   where
--     zippedToWindow = zip4 xs (drop 1 xs) (drop 2 xs) [1..]
--     keepValueAndIx = \(_, m, _, ix) -> (m, ix)
--     isLocalMaximum = \(l, m, r, _) -> m > l && m > r
--
--
-- findLocMinWithIx :: (Ord a, Num a) => [a] -> [(a, Int)]
-- findLocMinWithIx xs = map keepValueAndIx . filter isLocalMinimum $ zippedToWindow
--   where
--     zippedToWindow = zip4 xs (drop 1 xs) (drop 2 xs) [1..]
--     keepValueAndIx = \(_, m, _, ix) -> (m, ix)
--     isLocalMinimum = \(l, m, r, _) -> m < l && m < r


findLocMinMaxWithIx :: (Ord a, Num a) => [a] -> [(a, Int)]
findLocMinMaxWithIx xs = map keepValueAndIx . filter isLocalMinOrMax $ zippedToWindow
  where
    zippedToWindow  = zip4 xs (drop 1 xs) (drop 2 xs) [1..]
    keepValueAndIx  = \(_, m, _, ix) -> (m, ix)
    isLocalMinOrMax = \(l, m, r, _)  -> (m < l && m < r) || (m > l && m > r)


writePlot :: (Default r, ToRenderable r) => FilePath -> EC r () -> IO ()
writePlot filePath plot = toFile def {_fo_format=PDF} filePath plot
