{-# LANGUAGE NoImplicitPrelude    #-} -- disable standard prelude
{-# LANGUAGE RecordWildCards      #-} -- for passing Pars to ODEs
{-# LANGUAGE OverloadedStrings    #-} -- for Foundation
{-# LANGUAGE ExtendedDefaultRules #-} -- for Matplotlib

module Chemostat where

-- custom functions
import Util

-- instead of prelude
import Foundation
import Foundation.Collection ((!))

-- solving ODEs and working with the resulting data
import qualified Numeric.LinearAlgebra as LA
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)

-- plotting
import Graphics.Matplotlib

-- to keep type signatures short
type D = Double

-- parameter data type, holding the parameters that are passed to the differential equations.
data Par = Par { d, si, gc, xc, kmax, pu, pd, a, gp, kp, xp, gtp, ktp, xtp :: D }
  deriving (Show)

-- the parameters that are actually used in the simulation
basePars :: Par
basePars = Par {
  d    = 0.05, -- dilution rate
  si   = 4.0,  -- inflow nutrient concentration
  gc   = 0.5,  -- maximum consumption rate Cu/Cd
  xc   = 0.3,  -- conversion efficiency Cu/Cd 
  kmax = 0.3,  -- maximim half saturation Cu/Cd
  pu   = 1.0,  -- palatability Cu (undefended)
  pd   = 0.1,  -- palatability Cd (defended)
  a    = 0.29, -- trade-off strength
  gp   = 0.5,  -- maximum consumption rate Pi
  kp   = 0.5,  -- half saturation Pi
  xp   = 0.3,  -- conversion efficiency Pi
  gtp  = 0.5,  -- maximum consumption rate Pt
  ktp  = 0.05, -- half saturation Pt
  xtp  = 0.3   -- conversion effienciency Pt
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

-- maybe define the functional responses also as sub functions?

-- the differential equations to solve
ds, dcu, dcd, dp, dtp :: Par -> D -> D -> D -> D -> D -> D
ds  Par{..} s c1 c2 p tp = d*(si - s) - s*(gc*c1/(kcu + s) + (gc*c2/(kcd + s)))
dcu Par{..} s c1 c2 p tp = c1*(xc*gc*s/(kcu + s) - gp*pu*p/(kp + pu*c1 + pd*c2) - d)
dcd Par{..} s c1 c2 p tp = c2*(xc*gc*s/(kcd + s) - gp*pd*p/(kp + pu*c1 + pd*c2) - d)
dp  Par{..} s c1 c2 p tp = p*(xp*gp*((pu*c1 + pd*c2)/(kp + pd*c1 + pd*c2)) - gtp*tp/(ktp + p) - d)
dtp Par{..} s c1 c2 p tp = tp*(xtp*gtp*p/(ktp + p) - d)

-- the differential equations system in the form that is passed to the solver
-- pars need to be passed here for bifurcation to work later on
eqSystem :: Par -> D -> LA.Vector D -> LA.Vector D
eqSystem pars t vars = LA.fromList [ ds  basePars s c1 c2 p tp
                                   , dcu pars s c1 c2 p tp
                                   , dcd pars s c1 c2 p tp
                                   , dp  pars s c1 c2 p tp
                                   , dtp pars s c1 c2 p tp ]
  where 
    -- pattern matching on vars to get the population densities for using them in the functions body
    [s, c1, c2, p, tp] = LA.toList vars 

-- the time steps for which the result is given
time ::  LA.Vector D
time = times 0 10 1.0

-- solving the equations numerically; returning the solutions matrix
solveEqs :: (Par -> D -> LA.Vector D -> LA.Vector D) -> Par -> LA.Matrix D
solveEqs eqSys pars = odeSolveV
  RKf45    -- ODE Method
  1E-8     -- initial step size
  1E-8     -- absolute tolerance for the state vector
  0        -- relative tolerance for the state vector
  (eqSys pars) -- differential eqations: xdot(t,x), ...
  initVals -- inital conditions
  time     -- desired solution times

-- adding a column to the solution matrix, containing the total of both clones Cu+Cd
solWithCloneTotal :: LA.Matrix D
solWithCloneTotal = matrixAddSumCol 1 2 $ solveEqs eqSystem basePars

-- defining a plot that is later saved to a file
-- this uses python with matplotlib under the hood
plot1 :: Matplotlib
plot1 =
  plot time (LA.toColumns solWithCloneTotal ! 0) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "S"] %
  plot time (LA.toColumns solWithCloneTotal ! 1) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Cu"] %
  plot time (LA.toColumns solWithCloneTotal ! 2) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Cd"] %
  plot time (LA.toColumns solWithCloneTotal ! 3) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Pi"] %
  plot time (LA.toColumns solWithCloneTotal ! 4) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Pt"] %
  plot time (LA.toColumns solWithCloneTotal ! 5) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Ct"] %
  legend @@ [o2 "fancybox" True, o2 "shadow" False, o2 "title" "Legend", o2 "loc" "upper left"]

-- putting everyting together to export it for usage in Main.hs
runChemostat :: IO ()
runChemostat = do
  Right _ <- file "plot1.pdf" plot1
  return ()

-- start working on bifurcation, should maybe be refactored later
-- run time series, bifurcation and plots from different modules in Main.hs

bifurcationPars :: [Par]
bifurcationPars = [ basePars { d = x } | x <- [0.01,0.02..0.05] ]

-- bifEqSystemsWithPars :: [D -> LA.Vector D -> LA.Vector D]
-- bifEqSystemsWithPars = fmap eqSystem bifurcationPars

-- can i compute the systems and also output the currently processed par?
-- can i carry the used parameters with the computation, to later write some meta data?
bifSolutions :: [LA.Matrix D]
bifSolutions = fmap (solveEqs eqSystem) bifurcationPars

