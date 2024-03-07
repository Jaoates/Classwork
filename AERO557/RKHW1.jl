using Pkg
Pkg.activate("A557")
using Revise
includet("../orbitslib/OrbitsLibOD.jl")
includet("../orbitslib/OrbitsLibBase.jl")
using .OrbitsLibOD
using .OrbitsLibBase
using SatelliteToolbox
using LinearAlgebra
using Unzip
using Printf
using AstroLib
using GLMakie

##
## Prob 1

HST_tle = tle"""
HST                     
1 20580U 90037B   24010.28323428  .00005912  00000-0  29430-3 0  9991
2 20580  28.4707 327.8603 0002613 104.9653 010.4281 15.15476605652889
"""

lat, lon, h = [35.3540, -120.3757, 105.8]

dt_start = DateTime(2024, 1, 11) + Hour(8)  #UTC
dt_end = dt_start + Day(1)

jd = date_to_jd(dt_start):5/(3600*24):date_to_jd(dt_end)

orbp = Propagators.init(Val(:SGP4), HST_tle)
ret = Propagators.propagate_to_epoch!.(orbp, jd)
r = first.(ret)
r_ECEF = r_eci_to_ecef.(TEME(), PEF(), jd) .* r
r_NED = ecef_to_ned.(r_ECEF, lat |> deg2rad, lon |> deg2rad, h; translate=true)

α, δ = unzip(calcRADec.(r))
El, Az = unzip(calcElAz.(r_NED))

visible_idx = findall(El .> 0)
El_vis = El[visible_idx]
Az_vis = Az[visible_idx]
time_vis = jd_to_date.(DateTime, jd[visible_idx])

diff, indx = findmin(el -> abs(el - 10), El_vis)


# printstyled("Declination is $δ deg\n"; bold=true)
# printstyled("RA is $α deg\n"; bold=true)
# printstyled("Elevation is $El deg\n"; bold=true)
# printstyled("Azimuth is $Az deg\n"; bold=true)


##
## Problem 2
##
using Roots

RA = [30.381, 65.134, 99.976] #deg
DEC = [23.525, .774, -30.44]
date = "03-25-2013"
times = ["03:10:30.032", "03:15:20.612", "03:20:32.777"]
dts = [DateTime("$(date)T$time", dateformat"mm-dd-yyyyTHH:MM:SS.s") for time in times]
jd = date_to_jd.(dts)
lat, lon, h = [35.30, -120.66, 105.8] #deg deg m

r_site_ecef = geodetic_to_ecef(lat |> deg2rad, lon |> deg2rad, h)
r_site = r_ecef_to_eci.(PEF(), TEME(), jd) .* Ref(r_site_ecef)*1e-3 #km

ρ_h = [[cosd(δ)*cosd(α), cosd(δ)*sind(α), sind(δ)] for (α, δ) in zip(RA, DEC)]

v2_g, r2_g = gaussIOD(r_site, ρ_h, jd; ext=false)
kep_g = OrbitStateVector(jd[2], r2_g*1e3, v2_g*1e3) |> sv_to_kepler
v2_gext, r2_gext = gaussIOD(r_site, ρ_h, jd)
kep_gext = OrbitStateVector(jd[2], r2_gext*1e3, v2_gext*1e3) |> sv_to_kepler

using MATLAB
mat"addpath('mat')"
p1, p2, p3 = ρ_h
r1, r2, r3 = Vector.(r_site)
τ1 = (jd[1] - jd[2])*24*3600
τ3 = (jd[3] - jd[2])*24*3600
mat"""
[r2vec, v2vec] = AERO557doubleR($p1, $p2, $p3, $r1, $r2, $r3, $τ1, $τ3)
"""
r2_dr = @mget r2vec
v2_dr = @mget v2vec
kep_dr = OrbitStateVector(jd[2], r2_dr*1e3, v2_dr*1e3) |> sv_to_kepler


TLEs = tles"""
TLE1
1 25623U 99004C   13109.04882318 -.00000094  00000-0 -16690-3 0  1204
2 25623 051.9974 067.7982 0012092 184.1231 215.4516 11.86494224651645
TLE2
1 25165U 98008D   13083.14572197 -.00000211  00000-0 -12941-2 0  4434
2 25165 052.0160 303.6990 0005433 319.9573 182.6276 12.12023409691559
TLE3
1 25946U 99058D   13104.16396495 -.00000071  00000-0  19175-3 0  0939
2 25946 051.9981 329.5396 0000986 149.5293 353.4996 12.46940793622716
"""
for tle in TLEs
    orbp = Propagators.init(Val(:SGP4), tle)
    Propagators.propagate_to_epoch!(orbp, jd[2])
    kep = Propagators.mean_elements(orbp) 
    display(kep)
end

# using GLMakie
# # F = Figure()
# # azx = Axis3(F)
# Earth = Sphere(Point3f(0), rₑ)
# F, scene = mesh(Earth; alpha=.25)
# scatter!(scene, Point3f.(r_site); markersize=15)
# scatter!(scene, Point3f.(r); markersize=15) 


##
## Problem 3
##

r_0 = [15945.34, 0, 0]
r_1 = [12214.83899, 10249.46731, 0]
Δt = 76*60.0

# UV Method Lamberts 

v_UV_1, v_UV_2 = lamberts(r_0, r_1, Δt; tol=1e-12)
COEslam = calcCOE.([v_UV_1, v_UV_2], [r_0, r_1])
SV = OrbitStateVector.(Ref(0), 1e3*[r_0, r_1], 1e3*[v_UV_1, v_UV_2])
kepel_UV = sv_to_kepler.(SV)

# Izzo Method Lamberts
using PyCall
import Conda
pk = pyimport("pykep")
np = pyimport("numpy")

l = pk.lambert_problem(r1=r_0*1e3, r2=r_1*1e3, tof=Δt, mu=pk.MU_EARTH, max_revs=0)
v_izzo_1 = collect(l.get_v1()[1])*1e-3
v_izzo_2 = collect(l.get_v2()[1])*1e-3
kepel_izzo = sv_to_kepler.(OrbitStateVector.(Ref(0), 1e3*[r_0, r_1], 1e3*[v_izzo_1, v_izzo_2]))

# Gauss Method Lamberts


v_gauss_1 = lamberts_gauss(r_0, r_1, Δt)

# Minimum energy method lamberts 
μ = μₑ
r1 = r_0; r2 = r_1
Δθ = acosd(dot(r1, r2)/(norm(r1) * norm(r2)))
c = sqrt(norm(r1)^2 + norm(r2)^2 - 2*norm(r1)*norm(r2)*cosd(Δθ))
s = (norm(r1) + norm(r2) + c)/2
a_min = s/2
p_min = (1 - cosd(Δθ))*norm(r1)*norm(r2)/c
e_min = sqrt(1 - 2*p_min/s)
α_e = π
β_e = 2*asin(sqrt((s - c)/s))
t_min_e = sqrt(a_min^3/μ).*(α_e .+ (-1, 1).*(sin(β_e)))
t_min_abs = (1/3)*sqrt(2/μ) * (s^(1.5) - (s - c)^(1.5))
v_1 = (r2 .- (1 - norm(r2)/p_min*(1 - cosd(Δθ)))*r1 ) * sqrt(μ*p_min)/(norm(r1)*norm(r2)*sind(Δθ))