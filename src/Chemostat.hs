{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Chemostat (runChemostat) where

import Foundation
import Data.List ((!!))
-- import Numeric  ((**), exp)

import Graphics.Rendering.Chart.Easy ((.=), layout_title, plot, line)

import qualified Numeric.LinearAlgebra as LA
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)

import Util (getNowTimeString, makePlottableTuples, writePlot)

-- parameters
data Par = Par { d, si, gc, xc, kmax, p1, p2, a, gp, kp, xp, gtp, ktp, xtp :: Double }
  deriving (Show)

par :: Par
par = Par {
  d    = 0.05,
  si   = 4.0,
  gc   = 0.5,
  xc   = 0.3,
  kmax = 0.5,
  p1   = 1.0,
  p2   = 0.1,
  a    = 0.29,
  gp   = 0.5,
  kp   = 0.5,
  xp   = 0.3,
  gtp  = 0.5,
  ktp  = 0.05,
  xtp  = 0.3
}

-- initial values
initVals :: LA.Vector Double
initVals = LA.fromList [ 0.5, 0.5, 0.5, 0.2, 0.0 ]

-- sub functions
kc1', kc2' :: Par -> Double
kc1' Par{..} = kmax - a*p1
kc2' Par{..} = kmax - a*p2

kc1, kc2 :: Double
kc1 = kc1' par
kc2 = kc2' par

-- differential equations
ds, dc1, dc2, dp, dtp :: Par -> Double -> Double -> Double -> Double -> Double -> Double
ds   Par{..} s c1 c2 p tp = d*(si - s) - s*(gc*c1/(kc1 + s) + (gc*c2/(kc2 + s)))
dc1  Par{..} s c1 c2 p tp = c1*(xc*gc*s/(kc1 + s) - gp*p1*p/(kp + p1*c1 + p2*c2) - d)
dc2  Par{..} s c1 c2 p tp = c2*(xc*gc*s/(kc2 + s) - gp*p1*p/(kp + p1*c1 + p2*c2) - d)
dp   Par{..} s c1 c2 p tp = p*(xp*gp*((p1*c1 + p2*c2)/(kp + p1*c1 + p2*c2)) - gtp*tp/(ktp + p) - d)
dtp  Par{..} s c1 c2 p tp = tp*(xtp*gtp*p/(ktp + p) - d)

eqSystem :: Double -> LA.Vector Double -> LA.Vector Double
eqSystem t vars = LA.fromList [ ds  par s c1 c2 p tp
                              , dc1 par s c1 c2 p tp
                              , dc2 par s c1 c2 p tp
                              , dp  par s c1 c2 p tp
                              , dtp par s c1 c2 p tp ]
  where
    [s, c1, c2, p, tp] = LA.toList vars

-- the time steps for which the result is given
times :: LA.Vector Double
times = LA.linspace 1000 ( 0, 999 :: Double )

-- the solutions matrix
solution :: LA.Matrix Double
solution = odeSolveV
  RKf45    -- ODE Method
  1E-8     -- initial step size
  1E-8     -- absolute tolerance for the state vector
  0        -- relative tolerance for the state vector
  eqSystem -- differential eqations: xdot(t,x), ...
  initVals -- inital conditions [ x0, y0, α0, β0 ]
  times    -- desired solution times

timePlot = do
  layout_title .= "time series"
  plot $ line "substrate"    [makePlottableTuples times solution !! 0]
  plot $ line "clone 1"      [makePlottableTuples times solution !! 1]
  plot $ line "clone 2"      [makePlottableTuples times solution !! 2]
  plot $ line "predator"     [makePlottableTuples times solution !! 3]
  plot $ line "top-predator" [makePlottableTuples times solution !! 4]

-- phasePlot = do
  -- layout_title .= "Mougi / Iwasa – phase space"
  -- plot $ line "prey - predator" [ fmap (\ [x, y, _, _] -> (x, y)) $ LA.toLists solution ]

runChemostat :: IO ()
runChemostat = do
  timeStr <- getNowTimeString
  writePlot ("plots/chemo_timePlot_" <> timeStr <> ".pdf") timePlot
  -- writePlot ("plots/chemo_phasePlot_" <> timeStr <> ".pdf") phasePlot
