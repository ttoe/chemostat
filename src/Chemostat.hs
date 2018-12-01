module Chemostat where

import Control.Monad (forM)
import Data.List (nubBy)
import Flow
import Graphics.Rendering.Chart.Easy (EC, Layout, line, plot)
import Model0 (Par(..), basePars0, model0)
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)
import Numeric.LinearAlgebra
  ( Matrix
  , Vector
  , fromList
  , toColumns
  , toList
  , toLists
  )
import Util

type D = Double
type VecD = Vector Double
type MatD = Matrix Double

initVals :: VecD -- [  N,  Cu,  Cd,  Pg,  Ps]
initVals = fromList [0.5, 0.5, 0.5, 0.2, 0.0]

time :: VecD
time = times 0 10000 0.1

-- TODO: put in Model0; no need to export model itself, rather export model applied to solver
solveEqs :: (Par -> D -> VecD -> VecD) -> Par -> MatD
solveEqs model pars =
  odeSolveV
    RKf45 -- ODE Method
    1E-8 -- initial step size
    1E-8 -- absolute tolerance for the state vector
    1E-8 -- relative tolerance for the state vector
    (model pars) -- differential equations: xdot(t,x), ...
    initVals -- inital conditions
    time -- desired solution times

-- adding a column to the solution matrix, containing the total of both clones Cu+Cd, and the time
solWithTimeAndCt :: MatD
solWithTimeAndCt = matrixAddTimeCol time <| matrixAddSumCol 1 2 sol
  where
    sol = solveEqs model0 basePars0

-- TODO: switch to keymap
-- TODO: make a function of arguments from/to/step in Util
bifurcationPars :: [Par]
bifurcationPars = [basePars0 {d = x} | x <- [0.20,0.21 .. 0.30]]
  -- the parameters basically still exist in <bifurcationPars>
  -- with sequential evaluation they would be reusable

-- TODO: carry the used parameters with the computation, to later write some meta data?
-- bifurcation with progress ouput
bifSolutionsWO :: (Par -> D -> VecD -> VecD) -> [Par] -> IO [MatD]
bifSolutionsWO eqSystem bifPars = do
  let eqSolver = solveEqs eqSystem
  forM bifPars $ \pars -> do
    print pars
    pure <| eqSolver pars

-- bifurcation without progress output
bifSolutions :: [Par] -> [MatD]
bifSolutions bifPars = fmap (solveEqs model0) bifPars

timePlot sol = do
  plotLine "Cu" 2
  plotLine "Cd" 3
  plotLine "Pg" 4
  plotLine "Ps" 5
  plotLine "Ct" 6
  where
    plotLine name col = plot <| line name [mkPlotTuples sol !! col]

phasePlot1, phasePlot2 :: MatD -> EC (Layout D D) ()
phasePlot1 sol =
  plot <| line "Pg~Ct" [(\[_, _, _, _, pg, _, ct] -> (ct, pg)) <$> toLists sol]

phasePlot2 sol =
  plot <| line "Cu~Cd" [(\[_, cu, cd, _, _, _, _] -> (cu, cd)) <$> toLists sol]

runChemostat :: IO ()
runChemostat = do
  let addedSumColMatrices = matrixAddSumCol 1 2 <$> bifSolutions bifurcationPars
  let minsMaxs =
        fmap
          (fmap (findLocMinMax . toList) . toColumns . mtxLast 300)
          addedSumColMatrices
  let numMinsMaxs d = fmap (fmap $ length . nubBy (\x y -> x - y < d)) minsMaxs
  -- print minsMaxs
  print $ numMinsMaxs 1.0
  -- timeStr <- nowTimeString
  -- writePlot ("plots/chemo_timePlot_"  <> timeStr <> ".pdf") $ timePlot $ mtxLast 300 solWithTimeAndCt
  -- writePlot ("plots/chemo_phasePlot1_" <> timeStr <> ".pdf") $ phasePlot1 solWithTimeAndCt
  -- writePlot ("plots/chemo_phasePlot2_" <> timeStr <> ".pdf") $ phasePlot2 solWithTimeAndCt
