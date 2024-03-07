# PROBLEM 1
include("JoshOrbitsLib.jl")
using SatelliteToolbox
using Plots
using ReferenceFrameRotations
using DifferentialEquations

# actual probelm
δr₀ = [0,0,0]
δ̇r₀ = [1,-1,1]



T = 4 # units of time
n = 2π/T
tf = .25*T
_,δṙ= cwEquations(δr₀,δ̇r₀,tf,n)
# println(n*tf)
# println(δr)
# println(δṙ)

println("PROBLEM 1")
println("The final velocity of the particle relative to the origin at 1/4 of an orbit in m/s is:")
println(δṙ)

# PROBLEM 2
n = tdrift = NaN
# n = 1.002701 #revs/day
n = 1.00 #revs/day
n = n*2*pi #rad/day
n = n/(24*3600) #rads/s

δr₀ = [0,0,0] # intitial position km
δr = [-1,1,0]*10 # final drift position

tdrift = 2*3600
Φrr,Φrv,Φvr,Φvv = cwCoefs(tdrift,n)

# δr = Φrr * δr₀ + Φrv * δ̇r₀
# δṙ = Φvr * δr₀ + Φvv * δ̇r₀
δṙ₀⁺ = inv(Φrv)*(δr - Φrr * δr₀)
# δṙ⁻ = Φvr * δr₀ + Φvv * δ̇r₀⁺

δr⁻,δṙEndDrift = cwEquations(δr₀,δṙ₀⁺,tdrift,n)

tf = 6*3600 # 2hrs in sec
Φrr,Φrv,Φvr,Φvv = cwCoefs(tf,n)

δr₀ = [-1,1,0]*10 # intitial position km

δṙ⁺ = [0,0,0] # desired final velo
δr = [0,0,0] # desired final position

δṙ₀⁺ = inv(Φrv)*(δr - Φrr*δr₀)
δr⁻,δṙ⁻ = cwEquations(δr₀,δṙ₀⁺,tf,n)

Δv₀ = δṙ₀⁺ - δṙEndDrift
Δv = δṙ⁺ - δṙ⁻

ΔvTotal = (norm(Δv₀)+norm(Δv))*1e3

println("PROBLEM 2")
println("The Δv of the maneuver in m/s:")
println(ΔvTotal)


# Problem 3
z = 300
r = z+re
δr₀ = [-1,0,0]

_,n = calcPeriod(r)
tf = 30*60 # 30min
Φrr,Φrv,Φvr,Φvv = cwCoefs(tf,n)

δṙ⁺ = [0,0,0] # desired final velo
δr = [0,0,0] # desired final position

δṙ₀⁺ = inv(Φrv)*(δr - Φrr*δr₀)
δr⁻,δṙ⁻ = cwEquations(δr₀,δṙ₀⁺,tf,n)

Δv₀ = δṙ₀⁺ .- 0
Δv = δṙ⁺ - δṙ⁻

println("PROBLEM 3")
println("The Δv of the maneuver in m/s:")
ΔvTotal = (norm(Δv₀)+norm(Δv))*1e3
println(ΔvTotal)

# Problem 4
ra = rb = va = vb = NaN

ra = 8000
rb = 7000
va = calcSpeedCirc(ra)
vb = calcSpeedCirc(rb)

ra = [1,0,0]*ra
rb = [1,0,0]*rb

va = [0,1,0]*va
vb = [0,1,0]*vb

rRel,vRel,aRel =calcRelativeKinematics(ra,va,rb,vb)
speedRel = norm(vRel)

println("PROBLEM 4")
println("The relative speed in km/s is:")
println(speedRel)

