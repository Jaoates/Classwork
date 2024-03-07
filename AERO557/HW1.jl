include("JoshOrbitsLib2.jl")
# include("OtherOrbitsLib.jl")
using SatelliteToolbox
using Plots
using ReferenceFrameRotations
using DifferentialEquations
using Revise
using AstroLib

plotlyjs()

# Problem 1
# Basic path
# R ECI > Rsite > rho > AZ EL > RA DEC > is visible

# set up time
T0 = DateTime(2024,1,11,0)+ Hour(8)
Tf = DateTime(2024,1,11,24)+ Hour(8)
Tspan = [T0,Tf] .|> date_to_jd

npts = 10000 # number of points desired in orbit
Tspan = range(Tspan...,npts)

# set up TLE
tle = tle"""
1 20580U 90037B   24010.28323428  .00005912  00000-0  29430-3 0  9991
2 20580  28.4707 327.8603 0002613 104.9653 010.4281 15.15476605652889
"""

# build orbit
prop = Propagators.init(Val(:SGP4),tle)
elements = Propagators.mean_elements(prop)
# r,v = kepler_to_rv(elements)
# print(r)
# print(v)
prop = Propagators.init(Val(:TwoBody),elements)

# propagate orbits over tspan
rv = Propagators.propagate_to_epoch!.(prop,Tspan)
r = map(x->x[1]|>Vector,rv)
v = map(x->x[2]|>Vector,rv)
r/=1000
v/=1000

# Set up site location
Φ = 35.3540 |> deg2rad   # latitude
λ = -120.3757 |> deg2rad # longitude
h = 105.8                # altitude (m)

Rs = [Φ,λ,h]

# calc lst
# lst = ct2lst.(λ|>rad2deg,Tspan)|>deg2rad

# Rotate site location over time
Rs = geodetic_to_ecef(Rs...) # take R site to ECEF
q=r_ecef_to_eci.(Quaternion,PEF(),TEME(),Tspan,)
Rs = vecRotate.(Ref(Rs),q) # rotate into ECI frame

Rs ./= 1000 # km

# heart check
# print(maximum(abs,(re .- norm.(Rs))))
# plot(scrape.(Rs,1),scrape.(Rs,2),scrape.(Rs,3),show = true)
# plot!(aspect_ratio = :equal)
# axesEqual!()  

# rho
ρ = Rs .- r # rho in ECI

# 
αδ =  calcRaDec.(r)
α = scrape.(αδ,1)
δ = scrape.(αδ,2)

# elaz = calcElAz.(Φ,λ,α,Tspan)
elaz = map((x,y,z)->calcElAz(Φ,λ,x,y,z),α,δ,Tspan)

el = scrape.(elaz,1)
az = scrape.(elaz,2)

elCutoff = 0 |> deg2rad
elInd = el.>elCutoff
elInd2 = push!(diff(elInd),(diff(elInd)[end]))
elInd2 = elInd2
Tspan = Tspan |> Vector

Tstart = Tspan[elInd2 .== 1] 
Tend = Tspan[elInd2 .== -1]

# plot(Tspan,el|>Vector,show=true)

println("The satellite rises and sets with the horizon at the following times:")
println("number - start time - end time")
for i = 1:length(Tstart)
    println("$(i) - $(jd_to_date(DateTime,Tstart[i])) - $(jd_to_date(DateTime,Tend[i]))")
end

sun = sun_position_mod.(Tspan)./1000
Q = r_eci_to_eci.(Quaternion,MOD(),J2000(),Tspan)
sun = vecRotate.(sun,Q)

# find out if the satellite is in the sun
F = calcInShadow.(r,sun)
F = .!F
bothInd = elInd .& F
bothInd2 = push!(diff(bothInd),(diff(bothInd)[end]))

# find out if the site is in shadow
F2 = calcInShadow.(map(x->uvec(x).*(re+1),Rs),sun)
allInd = elInd .& F.& F2
allInd2 = push!(diff(allInd),(diff(allInd)[end]))

Tstart = Tspan[allInd2 .== 1]
Tend = Tspan[allInd2 .== -1]

println("The satellite is veiwable, and exits view at the following times:")
println("This means that the satellite is both above the horizon and in sun and the site is dark")
println("number - start time - end time")
for i = 1:length(Tstart)
    println("$(i) - $(jd_to_date(DateTime,Tstart[i])) - $(jd_to_date(DateTime,Tend[i]))")
end

scatter(Tspan,elInd)
plot!(Tspan,el,show = true)

# Problem 3
