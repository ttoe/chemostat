-- TODO: use (!) from foundation everywhere
-- TODO: use Debug.Trace (trace)

{-# LANGUAGE RecordWildCards      #-} -- for passing Pars to ODEs
{-# LANGUAGE OverloadedStrings    #-} -- for Foundation
{-# LANGUAGE ExtendedDefaultRules #-} -- for Matplotlib

module Chemostat where

-- timing
import Data.Time.Clock

-- custom functions
import Util

import Data.List ((!!))

-- instead of prelude
import Foundation
import Foundation.Collection ((!))

-- solving ODEs and working with the resulting data
import qualified Numeric.LinearAlgebra as LA
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)

import Control.Monad (forM)

-- plotting
import Graphics.Matplotlib

-- to keep type signatures short
type D = Double

-- parameter data type, holding the parameters that are passed to the differential equations.
data Par = Par { d, ni, gc, ec, kmax, pu, pd, a, gg, kg, eg, gs, ks, es :: D }
  deriving (Show)

basePars :: Par
basePars = Par {
  d    = 0.05, -- dilution rate
  ni   = 4.0,  -- inflow nutrient concentration
  gc   = 0.5,  -- maximum consumption rate Cu/Cd
  ec   = 0.3,  -- conversion efficiency Cu/Cd
  kmax = 0.3,  -- maximim half saturation Cu/Cd
  pu   = 1.0,  -- palatability Cu (undefended)
  pd   = 0.1,  -- palatability Cd (defended)
  a    = 0.29, -- trade-off strength
  gg   = 0.5,  -- maximum consumption rate Pi
  kg   = 0.5,  -- half saturation Pi
  eg   = 0.3,  -- conversion efficiency Pi
  gs   = 0.5,  -- maximum consumption rate Pt
  ks   = 0.05, -- half saturation Pt
  es   = 0.3   -- conversion effienciency Pt
}

-- initial values for all variables
initVals :: LA.Vector D -- [   S,  Cu,  Cd,  Pi,  Pt ]
initVals =  LA.fromList    [ 0.5, 0.5, 0.5, 0.2, 0.0 ]

-- helpers for sub-functions (half saturation/palatability trade-off)
-- or: kcu = let Par{..} = basePars in kmax - a * pu
-- or with ViewPatterns, PatternGuards, RecordWildcards even shorter:
-- kcu | Par{..} <- basePars = kmax - a*pu
kcu, kcd :: D
kcu = kmax - a*pu where Par{..} = basePars
kcd = kmax - a*pd where Par{..} = basePars

-- the differential equations to solve
dn, dcu, dcd, dpg, dps :: Par -> D -> D -> D -> D -> D -> D
dn  Par{..} n cu cd pg ps = - d*n
                            + d*ni
                            - n*gc*cu/(kcu + n)
                            - n*gc*cd/(kcd + n)
dcu Par{..} n cu cd pg ps = - d*cu
                            + cu*ec*gc*n/(kcu + n)
                            - cu*gg*pu*pg/(kg + pu*cu + pd*cd)
dcd Par{..} n cu cd pg ps = - d*cd
                            + cu*ec*gc*n/(kcd + n)
                            - gg*pd*pg/(kg + pu*cu + pd*cd)
dpg Par{..} n cu cd pg ps = - d*pg
                            + pg*eg*gg*(pu*cu)/(kg + pu*cu + pd*cd)
                            + pg*eg*gg*(pd*cd)/(kg + pu*cu + pd*cd)
                            - pg*gs*ps/(ks + pg)
dps Par{..} n cu cd pg ps = - d*ps
                            + ps*es*gs*pg/(ks + pg)

-- the differential equations system in the form that is passed to the solver
-- pars need to be passed here for bifurcation to work later on
eqSystem :: Par -> D -> LA.Vector D -> LA.Vector D
eqSystem pars t vars = LA.fromList [ dn  pars n cu cd pg ps
                                   , dcu pars n cu cd pg ps
                                   , dcd pars n cu cd pg ps
                                   , dpg pars n cu cd pg ps
                                   , dps pars n cu cd pg ps ]
  where
    -- pattern matching on vars to get the population densities for using them in the functions body
    [n, cu, cd, pg, ps] = LA.toList vars

-- the time steps for which the result is given
time ::  LA.Vector D
time = times 0 20000 0.1

-- solving the equations numerically; returning the solutions matrix
solveEqs :: (Par -> D -> LA.Vector D -> LA.Vector D) -> Par -> LA.Matrix D
solveEqs eqSys pars = odeSolveV
  RKf45        -- ODE Method
  1E-8         -- initial step size
  1E-8         -- absolute tolerance for the state vector
  0            -- relative tolerance for the state vector
  (eqSys pars) -- differential equations: xdot(t,x), ...
  initVals     -- inital conditions
  time         -- desired solution times

-- adding a column to the solution matrix, containing the total of both clones Cu+Cd
solWithCloneTotal :: LA.Matrix D
solWithCloneTotal = matrixAddSumCol 1 2 $ sol
  where
    sol = solveEqs eqSystem basePars

-- TODO: use frames instead of matrices and use semantic indexing rather than numbers

-- defining a plot that is later saved to a file
-- this uses python with matplotlib under the hood
myPlot :: Int -> LA.Vector D -> LA.Matrix D -> Matplotlib
myPlot drops time solMatrix =
  -- plot (cols ! 0) (cols ! 1) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "N"] %
  plot (cols !! 0) (cols !! 2) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Cu"] %
  plot (cols !! 0) (cols !! 3) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Cd"] %
  plot (cols !! 0) (cols !! 4) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Pg"] %
  plot (cols !! 0) (cols !! 5) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Ps"] %
  plot (cols !! 0) (cols !! 6) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Ct"] %
  legend @@ [o2 "fancybox" True, o2 "shadow" False, o2 "title" "Legend", o2 "loc" "upper left"]
  where
    withTime = LA.fromColumns $ time : LA.toColumns solMatrix
    sliced   = withTime LA.?? (LA.Drop drops, LA.All)
    cols     = LA.toColumns sliced

bifurcationPars :: [Par]
bifurcationPars = [ basePars { d = x } | x <- [0.01,0.0101..0.06] ]

-- bifEqSystemsWithPars :: [D -> LA.Vector D -> LA.Vector D]
-- bifEqSystemsWithPars = fmap eqSystem bifurcationPars

-- TODO: compute the systems and also output the currently processed par?
-- TODO: carry the used parameters with the computation, to later write some meta data?
  -- the parameters basically still exist in <bifurcationPars>
  -- with sequential evaluation they are reusable
-- Version without output:
  -- bifSolutions :: [Par] -> [LA.Matrix D]
  -- bifSolutions bifPars = fmap (solveEqs eqSystem) bifPars
bifSolutions :: [Par] -> IO [LA.Matrix D]
bifSolutions bifPars = do
  let eqSolver = solveEqs eqSystem

  res <- forM bifurcationPars $ \pars -> do
    putStrLn $ show (d pars)
    return $ eqSolver pars

  return res

addSumColToListOfMatrices :: Int -> Int -> [LA.Matrix D] -> [LA.Matrix D]
addSumColToListOfMatrices c1 c2 ms = fmap (matrixAddSumCol c1 c2) ms

-- putting everyting together to export it for usage in Main.hs
runChemostat :: IO ()
runChemostat = do
  startTime <- getCurrentTime
  -- sols <- mapM bifSolutions bifurcationPars
  sols <- bifSolutions bifurcationPars
  let addedColMatrices = addSumColToListOfMatrices 1 2 sols
  -- let addedColMatrices = fmap (matrixAddSumCol 1 2) sols
  endTime   <- getCurrentTime
  putStrLn $ show $ diffUTCTime endTime startTime
  -- Right _ <- file "plot0.pdf" $ myPlot solWithCloneTotal
  -- Right _ <- file "plot1.pdf" $ myPlot 10000 time $ addedColMatrices !! 0
  -- Right _ <- file "plot2.pdf" $ myPlot 10000 time $ addedColMatrices !! 100
  -- Right _ <- file "plot3.pdf" $ myPlot 10000 time $ addedColMatrices !! 200
  -- Right _ <- file "plot4.pdf" $ myPlot 10000 time $ addedColMatrices !! 300
  Right _ <- file "plot5.pdf" $ myPlot 197000 time $ addedColMatrices !! 500
  return ()
