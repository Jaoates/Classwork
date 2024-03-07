include("JoshOrbitsLib2.jl")
using SatelliteToolbox
using AstroLib

lat = 35.3540  |> deg2rad
lon = -120.3757 |>  deg2rad
jd = 2000000
Ra = π/4
Dec = π/4


lst = calcLst(lon,jd)
lst-lon
calcLocHourAngle(lst,Ra)

# calcElAz(lat,lon,Ra,Dec,jd)


