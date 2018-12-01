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
  ) where

import Data.List (zip3, zip4)
import Data.Time.Clock (UTCTime, diffUTCTime, getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)
import qualified Data.Vector.Storable as VS
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
type MatD = Matrix Double
type VecD = Vector Double

-- alias differentiation function to central derivative with initial step size 0.01
diff :: (D -> D) -> D -> D
diff fun point = fst $ derivCentral 0.01 fun point

-- takes a matrix where column 0 is the time
-- and the other columns are the solutions at a given time
-- creates a list of lists containing tuples
-- for each variable and the corresponding time
mkPlotTuples :: MatD -> [[(D, D)]]
mkPlotTuples m = fmap (zip timeCol) columns
  where
    columns = toList <$> toColumns m
    timeCol = head columns

mtxLast :: Int -> MatD -> MatD
mtxLast l m = m ?? (TakeLast takeLast, All)
  where
    timeStep = m `atIndex` (1, 0) -- 2nd value (1, _) in 1st column (_, 0) is time step from times function
    takeLast = floor $ fromIntegral l / timeStep

-- create a timestamp like "2017-06-09_131211"
nowTimeString :: IO String
nowTimeString = formatTime defaultTimeLocale "%F_%H%M%S" <$> getCurrentTime

matrixAddSumCol :: Int -> Int -> MatD -> MatD
matrixAddSumCol x y m = fromColumns (columns <> [colSum])
  where
    columns = toColumns m
    colSum = VS.zipWith (+) (columns !! x) (columns !! y)

matrixAddTimeCol :: VecD -> MatD -> MatD
matrixAddTimeCol t m = fromColumns ([t] <> toColumns m)

times :: D -> D -> D -> VecD
times from to stepsize = linspace steps (from, to)
  where
    steps = floor $ (to - from) / stepsize

printTimeDiff :: String -> UTCTime -> IO ()
printTimeDiff s startTime = do
  endTime <- getCurrentTime
  print $ s <> show (diffUTCTime endTime startTime)

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

findLocMinMax :: (Ord a, Num a) => [a] -> [a]
findLocMinMax xs = map keepValue . filter isLocalMinOrMax $ zippedToWindow
  where
    zippedToWindow = zip3 xs (drop 1 xs) (drop 2 xs)
    keepValue (_, m, _) = m
    isLocalMinOrMax (l, m, r) = (m < l && m < r) || (m > l && m > r)

writePlot :: (Default r, ToRenderable r) => FilePath -> EC r () -> IO ()
writePlot filePath plot = toFile def {_fo_format = PDF} filePath plot
