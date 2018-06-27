{-# LANGUAGE RecordWildCards #-} -- passing data type fields to functions

module Chemostat where

import Numeric.LinearAlgebra (Matrix, Vector, fromList, toList, toLists)
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)
import Control.Monad (forM)
import qualified Graphics.Rendering.Chart.Easy as P
import Util


type D = Double


-- parameter data type, holding the parameters that are passed to the differential equations.
data Par = Par { d, ni, r, pu, pd, ku, kd, gg, kg, eg, gs, ks, es :: D }
  deriving (Show)


basePars :: Par
basePars = Par {
  d    = 0.25, -- dilution rate
  ni   = 80.0, -- inflow nutrient concentration
  r    = 1.0,  -- growth rate clones
  pu   = 1.0,  -- palatability Cu (undefended)
  pd   = 0.01, -- palatability Cd (defended)
  ku   = 4.3,  -- half saturation Cu
  kd   = 8.6,  -- half saturation Cd
  gg   = 1.76, -- maximum consumption rate Pg
  kg   = 28.0, -- half saturation Pg
  eg   = 0.5,  -- conversion efficiency Pg
  gs   = 1.76, -- maximum consumption rate Ps
  ks   = 2.0,  -- half saturation Ps
  es   = 0.5   -- conversion effienciency Ps
}


-- initial values for all variables
initVals :: Vector D -- [   N,  Cu,  Cd,  Pg,  Ps ]
initVals =  fromList    [ 0.5, 0.5, 0.5, 0.2, 0.0 ]


{-
-- helpers for sub-functions (half saturation/palatability trade-off)
kcu, kcd :: D
kcu = kmax - a*pu where Par{..} = basePars
kcd = kmax - a*pd where Par{..} = basePars
-- kcu = let Par{..} = basePars in kmax - a * pu
-}


-- the differential equations to solve
dn, dcu, dcd, dpg, dps :: Par -> D -> D -> D -> D -> D -> D
dn  Par{..} n cu cd _  _  = - d*n
                            + d*ni
                            - r*(n/(ku+n))*cu
                            - r*(n/(kd+n))*cd
dcu Par{..} n cu cd pg _  = - d*cu
                            + r*(n/(ku+n))*cu
                            - gg*pu*(cu/(kg + pu*cu + pd*cd))*pg
dcd Par{..} n cu cd pg _  = - d*cd
                            + r*(n/(kd+n))*cd
                            - gg*pd*(cd/(kg + pu*cu + pd*cd))*pg
dpg Par{..} _ cu cd pg ps = - d*pg
                            + eg*gg*pu*(cu/(kg + pu*cu + pd*cd))*pg
                            + eg*gg*pd*(cd/(kg + pu*cu + pd*cd))*pg
                            - gs*(pg/(ks+pg))*ps
dps Par{..} _ _  _  pg ps = - d*ps
                            + es*gs*(pg/(ks + pg))*ps


-- the differential equations system in the form that is passed to the solver
eqSystem :: Par -> D -> Vector D -> Vector D
eqSystem pars t vars = fromList [ dn  pars n cu cd pg ps
                                , dcu pars n cu cd pg ps
                                , dcd pars n cu cd pg ps
                                , dpg pars n cu cd pg ps
                                , dps pars n cu cd pg ps ]
  where
    -- pattern matching on vars to get the population densities for using them in the functions body
    [n, cu, cd, pg, ps] = toList vars


-- the time steps for which the result is given
time :: Vector D
time = times 0 1000 0.1


-- solving the equations numerically; returning the solutions matrix
solveEqs :: (Par -> D -> Vector D -> Vector D) -> Par -> Matrix D
solveEqs eqSys pars = odeSolveV
  RKf45        -- ODE Method
  1E-8         -- initial step size
  1E-8         -- absolute tolerance for the state vector
  0            -- relative tolerance for the state vector
  (eqSys pars) -- differential equations: xdot(t,x), ...
  initVals     -- inital conditions
  time         -- desired solution times


-- adding a column to the solution matrix, containing the total of both clones Cu+Cd
solWithCt :: Matrix D
solWithCt = matrixAddSumCol 1 2 $ sol
  where
    sol = solveEqs eqSystem basePars


bifurcationPars :: [Par]
bifurcationPars = [ basePars { d = x } | x <- [0.01,0.015..0.06] ]


-- TODO: carry the used parameters with the computation, to later write some meta data?
  -- the parameters basically still exist in <bifurcationPars>
  -- with sequential evaluation they would be reusable
-- bifurcation with progress ouput
bifSolutionsWO :: [Par] -> IO [Matrix D]
bifSolutionsWO bifPars = do
  let eqSolver = solveEqs eqSystem

  res <- forM bifPars $ \pars -> do
    putStrLn $ show pars
    pure $ eqSolver pars
  pure res


-- bifurcation without progress output
bifSolutions :: [Par] -> [Matrix D]
bifSolutions bifPars = fmap (solveEqs eqSystem) bifPars


timePlot sol = do
  -- P.plot $ P.line "N"  [mkPlotTuples time sol !! 0]
  P.plot $ P.line "Cu" [mkPlotTuples time sol !! 1]
  P.plot $ P.line "Cd" [mkPlotTuples time sol !! 2]
  P.plot $ P.line "Pg" [mkPlotTuples time sol !! 3]
  P.plot $ P.line "Ps" [mkPlotTuples time sol !! 4]
  P.plot $ P.line "Ct" [mkPlotTuples time sol !! 5]


phasePlot sol = do
  P.plot $ P.line "Pg~Ct" [ fmap (\[_, _, _, pg, _, ct] -> (ct, pg)) $ toLists sol ]


runChemostat :: IO ()
runChemostat = do
  let addedSumColMatrices = fmap (matrixAddSumCol 1 2) $ bifSolutions bifurcationPars

  timeStr <- nowTimeString
  writePlot ("plots/chemo_timePlot_"  <> timeStr <> ".pdf") $ timePlot  solWithCt
  writePlot ("plots/chemo_phasePlot_" <> timeStr <> ".pdf") $ phasePlot solWithCt
