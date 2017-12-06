{-# LANGUAGE RecordWildCards #-} -- for passing Pars to ODEs
{-# LANGUAGE OverloadedStrings #-} -- for Foundation
{-# LANGUAGE NoImplicitPrelude #-} -- disable standard prelude
{-# LANGUAGE ExtendedDefaultRules #-} -- for Matplotlib

module Chemostat (runChemostat) where

import Util
import Foundation
import Foundation.Collection ((!))
import qualified Numeric.LinearAlgebra as LA
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)
import Graphics.Matplotlib

type D = Double

-- parameters
data Par = Par { d, si, gc, xc, kmax, pu, pd, a, gp, kp, xp, gtp, ktp, xtp :: D }

par :: Par
par = Par {
  d    = 0.05,
  si   = 4.0,
  gc   = 0.5,
  xc   = 0.3,
  kmax = 0.3,
  pu   = 1.0,
  pd   = 0.1,
  a    = 0.29,
  gp   = 0.5,
  kp   = 0.5,
  xp   = 0.3,
  gtp  = 0.5,
  ktp  = 0.05,
  xtp  = 0.3
} 

-- initial values
initVals :: LA.Vector D
initVals = LA.fromList [ 0.5, 0.5, 0.5, 0.2, 0.0 ]

-- sub functions
kcu', kcd' :: Par -> D
kcu' Par{..} = kmax - a*pu
kcd' Par{..} = kmax - a*pd

kcu, kcd :: D
kcu = kcu' par
kcd = kcd' par

-- differential equations
ds, dcu, dcd, dp, dtp :: Par -> D -> D -> D -> D -> D -> D
ds   Par{..} s c1 c2 p tp = d*(si - s) - s*(gc*c1/(kcu + s) + (gc*c2/(kcd + s)))
dcu  Par{..} s c1 c2 p tp = c1*(xc*gc*s/(kcu + s) - gp*pu*p/(kp + pu*c1 + pd*c2) - d)
dcd  Par{..} s c1 c2 p tp = c2*(xc*gc*s/(kcd + s) - gp*pd*p/(kp + pu*c1 + pd*c2) - d)
dp   Par{..} s c1 c2 p tp = p*(xp*gp*((pu*c1 + pd*c2)/(kp + pd*c1 + pd*c2)) - gtp*tp/(ktp + p) - d)
dtp  Par{..} s c1 c2 p tp = tp*(xtp*gtp*p/(ktp + p) - d)

eqSystem :: D -> LA.Vector D -> LA.Vector D
eqSystem t vars = LA.fromList [ ds  par s c1 c2 p tp
                              , dcu par s c1 c2 p tp
                              , dcd par s c1 c2 p tp
                              , dp  par s c1 c2 p tp
                              , dtp par s c1 c2 p tp ]
  where
    [s, c1, c2, p, tp] = LA.toList vars

-- the time steps for which the result is given
time ::  LA.Vector D
time = times 0 1000 0.1

-- the solutions matrix
solution :: LA.Matrix D
solution = odeSolveV
  RKf45    -- ODE Method
  1E-8     -- initial step size
  1E-8     -- absolute tolerance for the state vector
  0        -- relative tolerance for the state vector
  eqSystem -- differential eqations: xdot(t,x), ...
  initVals -- inital conditions
  time     -- desired solution times

solWithCloneTotal :: LA.Matrix D
solWithCloneTotal = matrixAddSumCol 1 2 solution

plot1 :: Matplotlib
plot1 =
  plot time (LA.toColumns solWithCloneTotal ! 0) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "S"] %
  plot time (LA.toColumns solWithCloneTotal ! 1) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Cu"] %
  plot time (LA.toColumns solWithCloneTotal ! 2) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Cd"] %
  plot time (LA.toColumns solWithCloneTotal ! 3) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Pi"] %
  plot time (LA.toColumns solWithCloneTotal ! 4) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Pt"] %
  plot time (LA.toColumns solWithCloneTotal ! 5) @@ [o1 "-", o2 "linewidth" 1, o2 "label" "Ct"] %
  legend @@ [o2 "fancybox" True, o2 "shadow" False, o2 "title" "Legend", o2 "loc" "upper left"]

runChemostat :: IO ()
runChemostat = do
  Right _ <- file "plot1.pdf" plot1
  return ()
