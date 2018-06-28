module Chemostat where

import Control.Monad (forM)
import Data.List (nubBy)
import qualified Graphics.Rendering.Chart.Easy as P
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)
import Numeric.LinearAlgebra
  ( Matrix
  , Vector
  , fromList
  , toColumns
  , toList
  , toLists
  )

import qualified Model0 as M0
import Util

type D = Double

initVals :: Vector D -- [   N,  Cu,  Cd,  Pg,  Ps ]
initVals = fromList [0.5, 0.5, 0.5, 0.2, 0.0]

time :: Vector D
time = times 0 10000 0.1

-- TODO: put in Model0; no need to export model itself, rather export model applied to solver
solveEqs :: (M0.Par -> D -> Vector D -> Vector D) -> M0.Par -> Matrix D
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
solWithTimeAndCt :: Matrix D
solWithTimeAndCt = matrixAddTimeCol time $ matrixAddSumCol 1 2 sol
  where
    sol = solveEqs M0.model M0.basePars

-- TODO: switch to keymap
-- TODO: make a function of arguments from/to/step in Util
bifurcationPars :: [M0.Par]
bifurcationPars = [M0.basePars {M0.d = x} | x <- [0.20,0.21 .. 0.30]]
  -- the parameters basically still exist in <bifurcationPars>
  -- with sequential evaluation they would be reusable

-- TODO: carry the used parameters with the computation, to later write some meta data?
-- bifurcation with progress ouput
bifSolutionsWO ::
     (M0.Par -> D -> Vector D -> Vector D) -> [M0.Par] -> IO [Matrix D]
bifSolutionsWO eqSystem bifPars = do
  let eqSolver = solveEqs eqSystem
  forM bifPars $ \pars -> do
    print pars
    pure $ eqSolver pars

-- bifurcation without progress output
bifSolutions :: [M0.Par] -> [Matrix D]
bifSolutions bifPars = fmap (solveEqs M0.model) bifPars

timePlot sol = do
  plotLine "Cu" 2
  plotLine "Cd" 3
  plotLine "Pg" 4
  plotLine "Ps" 5
  plotLine "Ct" 6
  where
    plotLine name col = P.plot $ P.line name [mkPlotTuples sol !! col]

phasePlot1, phasePlot2 :: Matrix D -> P.EC (P.Layout D D) ()
phasePlot1 sol =
  P.plot $
  P.line "Pg~Ct" [(\[_, _, _, _, pg, _, ct] -> (ct, pg)) <$> toLists sol]

phasePlot2 sol =
  P.plot $
  P.line "Cu~Cd" [(\[_, cu, cd, _, _, _, _] -> (cu, cd)) <$> toLists sol]

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
