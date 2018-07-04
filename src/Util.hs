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
  , findLocMinMax
  , findLocMinMaxV
  , flmm
  ) where

import Data.List (zip3, zip4)
import Data.Time.Clock (UTCTime, diffUTCTime, getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)
import qualified Data.Vector as V
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Storable as VS
import Flow
import Graphics.Rendering.Chart.Backend.Cairo
  ( FileFormat(..)
  , _fo_format
  , toFile
  )
import Graphics.Rendering.Chart.Easy (Default, EC, ToRenderable, def)
import Numeric.GSL.Differentiation (derivCentral)
import Numeric.LinearAlgebra
  ( Extractor(..)
  , Matrix
  , Vector
  , (??)
  , atIndex
  , fromColumns
  , linspace
  , toColumns
  , toList
  )

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
    columns = toList <$> toColumns m
    timeCol = head columns

mtxLast :: Int -> Matrix D -> Matrix D
mtxLast l m = m ?? (TakeLast takeLast, All)
  where
    timeStep = m `atIndex` (1, 0) -- 2nd value (1, _) in 1st column (_, 0) is time step from times function
    takeLast = floor $ fromIntegral l / timeStep

-- create a timestamp like "2017-06-09_131211"
nowTimeString :: IO String
nowTimeString = formatTime defaultTimeLocale "%F_%H%M%S" <$> getCurrentTime

matrixAddSumCol :: Int -> Int -> Matrix D -> Matrix D
matrixAddSumCol x y m = fromColumns (columns <> [colSum])
  where
    columns = toColumns m
    colSum = VS.zipWith (+) (columns !! x) (columns !! y)

matrixAddTimeCol :: Vector D -> Matrix D -> Matrix D
matrixAddTimeCol t m = fromColumns ([t] <> toColumns m)

times :: D -> D -> D -> Vector D
times from to stepsize = linspace steps (from, to)
  where
    steps = floor $ (to - from) / stepsize

printTimeDiff :: String -> UTCTime -> IO ()
printTimeDiff s startTime = do
  endTime <- getCurrentTime
  print $ s <> " | " <> show (diffUTCTime endTime startTime)

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
-- TODO: handle plateaus
findLocMinMaxWithIx :: (Ord a, Num a) => [a] -> [(a, Int)]
findLocMinMaxWithIx xs =
  map keepValueAndIx . filter isLocalMinOrMax $ zippedToWindow
  where
    zippedToWindow = zip4 xs (drop 1 xs) (drop 2 xs) [1 ..]
    keepValueAndIx (_, m, _, ix) = (m, ix)
    isLocalMinOrMax (l, m, r, _) = (m < l && m < r) || (m > l && m > r)

findLocMinMaxV :: VS.Vector D -> V.Vector D
findLocMinMaxV vec =
  V.map keepValue <| V.filter isLocalMinOrMax (zipped convVec)
  where
    convVec = V.convert vec :: V.Vector D
    zipped :: V.Vector D -> V.Vector (D, D, D)
    zipped v = V.zip3 v (V.drop 1 v) (V.drop 2 v)
    keepValue (_, m, _) = m
    isLocalMinOrMax (l, m, r) = (m < l && m < r) || (m > l && m > r)

findLocMinMax :: Ord a => [a] -> [a]
findLocMinMax xs = map keepValue . filter isLocalMinOrMax $ zippedToWindow
  where
    zippedToWindow = zip3 xs (drop 1 xs) (drop 2 xs)
    keepValue (_, m, _) = m
    isLocalMinOrMax (l, m, r) = (m < l && m < r) || (m > l && m > r)

flmm :: VS.Vector D -> VS.Vector D
flmm vec = VS.ifilter (\ix val -> boolVec vec VS.! ix) <| VS.init (VS.tail vec)
  where
    isLocalMinOrMax :: D -> D -> D -> Bool
    isLocalMinOrMax l m r = (m < l && m < r) || (m > l && m > r)
    boolVec :: VS.Vector D -> VS.Vector Bool
    boolVec v = VS.zipWith3 isLocalMinOrMax v (VS.drop 1 v) (VS.drop 2 v)

--   where
--     trueIx = ix + 1
--     lIx =
--     mIx =
--     rIx =
--     isLocalMinOrMax :: Int -> Bool
--     isLocalMinOrMax ix = (m < l && m < r) || (m > l && m > r)
writePlot :: (Default r, ToRenderable r) => FilePath -> EC r () -> IO ()
writePlot filePath plot = toFile def {_fo_format = PDF} filePath plot
