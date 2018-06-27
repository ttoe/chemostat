module Chemostat where

import Data.Time.Clock (getCurrentTime)
import Numeric.LinearAlgebra (Matrix, Vector, fromList, toList, toLists, toColumns)
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)
import Control.Monad (forM)
import qualified Graphics.Rendering.Chart.Easy as P

import Util
import qualified Model0 as M0


type D = Double


initVals :: Vector D -- [   N,  Cu,  Cd,  Pg,  Ps ]
initVals =  fromList    [ 0.5, 0.5, 0.5, 0.2, 0.0 ]


time :: Vector D
time = times 0 10000 0.1


-- TODO: put in Model0; no need to export model itself, rather export model applied to solver
solveEqs :: (M0.Par -> D -> Vector D -> Vector D) -> M0.Par -> Matrix D
solveEqs model pars = odeSolveV
  RKf45        -- ODE Method
  1E-8         -- initial step size
  1E-8         -- absolute tolerance for the state vector
  1E-8         -- relative tolerance for the state vector
  (model pars) -- differential equations: xdot(t,x), ...
  initVals     -- inital conditions
  time         -- desired solution times


-- adding a column to the solution matrix, containing the total of both clones Cu+Cd, and the time
solWithTimeAndCt :: Matrix D
solWithTimeAndCt = matrixAddTimeCol time $ matrixAddSumCol 1 2 $ sol
  where
    sol = solveEqs M0.model M0.basePars


-- TODO: switch to keymap
-- TODO: make a function of arguments from/to/step in Util
bifurcationPars :: [M0.Par]
bifurcationPars = [ M0.basePars { M0.d = x } | x <- [0.01,0.015..0.06] ]


-- TODO: carry the used parameters with the computation, to later write some meta data?
  -- the parameters basically still exist in <bifurcationPars>
  -- with sequential evaluation they would be reusable
-- bifurcation with progress ouput
bifSolutionsWO :: (M0.Par -> D -> Vector D -> Vector D) -> [M0.Par] -> IO [Matrix D]
bifSolutionsWO eqSystem bifPars = do
  let eqSolver = solveEqs eqSystem

  res <- forM bifPars $ \pars -> do
    putStrLn $ show pars
    pure $ eqSolver pars
  pure res


-- bifurcation without progress output
bifSolutions :: [M0.Par] -> [Matrix D]
bifSolutions bifPars = fmap (solveEqs M0.model) bifPars


timePlot sol = do
  -- plotLine "N" 1
  plotLine "Cu" 2
  plotLine "Cd" 3
  plotLine "Pg" 4
  plotLine "Ps" 5
  plotLine "Ct" 6
  where
    plotLine name col = P.plot $ P.line name [mkPlotTuples sol !! col]


phasePlot1 sol = do
  P.plot $ P.line "Pg~Ct" [ fmap (\[_, _, _, _, pg, _, ct] -> (ct, pg)) $ toLists sol ]

phasePlot2 sol = do
  P.plot $ P.line "Cu~Cd" [ fmap (\[_, cu, cd, _, _, _, _] -> (cu, cd)) $ toLists sol ]

runChemostat :: IO ()
runChemostat = do
  let addedSumColMatrices = fmap (matrixAddSumCol 1 2) $ bifSolutions bifurcationPars

  -- timeStr <- nowTimeString
  -- writePlot ("plots/chemo_timePlot_"  <> timeStr <> ".pdf") $ timePlot $ mtxLast 300 solWithTimeAndCt
  -- writePlot ("plots/chemo_phasePlot1_" <> timeStr <> ".pdf") $ phasePlot1 solWithTimeAndCt
  -- writePlot ("plots/chemo_phasePlot2_" <> timeStr <> ".pdf") $ phasePlot2 solWithTimeAndCt
