module Util
  ( diff
  , makePlottableTuples
  , getNowTimeString
  , writePlot
  ) where

import Foundation
import qualified Prelude as P

import qualified Graphics.Rendering.Chart.Backend.Cairo as BC
import Graphics.Rendering.Chart.Easy (def)

import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import qualified Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra (Vector, Matrix, toColumns)
import Numeric.GSL.Differentiation (derivCentral)

-- alias differentiation function to central derivative with initial step size 0.01
diff :: (Double -> Double) -> Double -> Double
diff fun point = fst $ derivCentral 0.01 fun point

-- takes the time vector and the solution matrix
-- creates a list of lists containing tuples
-- for each variable and the corresponding time
-- necessary to satisfy plotting funtion
makePlottableTuples :: Vector Double -> Matrix Double -> [[(Double, Double)]]
makePlottableTuples v m = fmap (P.zip (LA.toList v) . LA.toList) (toColumns m)

-- create a string like "2017-06-09_131211"
-- as a timestamp to use for writing files
getNowTimeString :: IO P.String
getNowTimeString = do
  now <- getCurrentTime
  pure (formatTime defaultTimeLocale "%F_%H%M%S" now)

writePlot filePath plot = BC.toFile def {BC._fo_format=BC.PDF} filePath plot
