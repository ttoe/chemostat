-- TODO: reexport Numeric.LinearAlgebra functions
{-# LANGUAGE RecordWildCards #-} -- passing data type fields to functions

module Model0 where

import Numeric.LinearAlgebra (Vector, fromList, toList)

type D = Double

-- parameter data type, holding the parameters that are passed to the differential equations.
data Par = Par
  { d, ni, r, pu, pd, ku, kd, gg, kg, eg, gs, ks, es :: D
  } deriving (Show)

basePars0 :: Par
basePars0 =
  Par
    { d = 0.25 -- dilution rate
    , ni = 80.0 -- inflow nutrient concentration
    , r = 1.0 -- growth rate clones
    , pu = 1.0 -- palatability Cu (undefended)
    , pd = 0.01 -- palatability Cd (defended)
    , ku = 4.3 -- half saturation Cu
    , kd = 8.6 -- half saturation Cd
    , gg = 1.76 -- maximum consumption rate Pg
    , kg = 28.0 -- half saturation Pg
    , eg = 0.5 -- conversion efficiency Pg
    , gs = 1.76 -- maximum consumption rate Ps
    , ks = 2.0 -- half saturation Ps
    , es = 0.5 -- conversion effienciency Ps
    }

dn, dcu, dcd, dpg, dps :: Par -> D -> D -> D -> D -> D -> D
dn Par {..} n cu cd _ _ =
  -d * n + d * ni - r * (n / (ku + n)) * cu - r * (n / (kd + n)) * cd

dcu Par {..} n cu cd pg _ =
  -d * cu + r * (n / (ku + n)) * cu -
  gg * pu * (cu / (kg + pu * cu + pd * cd)) * pg

dcd Par {..} n cu cd pg _ =
  -d * cd + r * (n / (kd + n)) * cd -
  gg * pd * (cd / (kg + pu * cu + pd * cd)) * pg

dpg Par {..} _ cu cd pg ps =
  -d * pg + eg * gg * pu * (cu / (kg + pu * cu + pd * cd)) * pg +
  eg * gg * pd * (cd / (kg + pu * cu + pd * cd)) * pg -
  gs * (pg / (ks + pg)) * ps

dps Par {..} _ _ _ pg ps = -d * ps + es * gs * (pg / (ks + pg)) * ps

model0 :: Par -> D -> Vector D -> Vector D
model0 pars t vars =
  fromList
    [ dn pars n cu cd pg ps
    , dcu pars n cu cd pg ps
    , dcd pars n cu cd pg ps
    , dpg pars n cu cd pg ps
    , dps pars n cu cd pg ps
    ]
    -- pattern matching on vars to get the population densities for using them in the functions body
  where
    [n, cu, cd, pg, ps] = toList vars
