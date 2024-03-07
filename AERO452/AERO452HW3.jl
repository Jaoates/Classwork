
include("JoshOrbitsLib.jl")
include("OtherOrbitsLib.jl")
using SatelliteToolbox
using Plots
# using ReferenceFrameRotations
using DifferentialEquations
using Revise

plotlyjs()

# Problem 1
println("--Problem 1--")
# Cowells
println("-Cowells-")
function EOM(X,p,t)
    μ = p.μ
    BC = p.BC
    r = X[1:3]
    v = X[4:6]
    a = calcAccel(r,μ) + calcDrag(r,v,BC)
    Ẋ = [v;a]
end
struct Params
    μ
    BC
    jd
end
Params(μ,BC)=Params(μ,BC,NaN)

m = 100
Cd = 2.2
A = π*(1/2)^2
BC = m/(Cd*A)
p = Params(μe,BC)

rp = 215 + re
ra = 939 + re
ecc = (ra-rp)/(ra+rp)
a = (ra+rp)/2

inc = 65.2 |> deg2rad
Ω = 340|> deg2rad
ω = 58|> deg2rad
θ = 332 |> deg2rad

T0 = date_to_jd(-4713, 11, 24, 12, 0, 0.0)

orb₀ = KeplerianElements(T0,a*1000,ecc,inc,Ω,ω,θ)
r₀,v₀ = kepler_to_rv(orb₀)
r₀ /=1000 
v₀ /=1000 

X₀ = [r₀;v₀] |> Vector

function myCondition(X,t,integrator)
    r = X[1:3]
    # v = X[4:6]
    (norm(r)-re)<100
end

myCB = DiscreteCallback(myCondition,terminate!)

tspan = (0,1000*24*3600)
myProb = ODEProblem(EOM,X₀,tspan,p)

sol = @time solve(myProb,reltol = 1e-8,abstol = 1e-8,callback = myCB)

r = [sol[1,:] sol[2,:] sol[3,:]]
v = [sol[4,:] sol[5,:] sol[6,:]]
t = sol.t
t /= (24*3600)
tf = t[end]
# plot(eachcol(r)...,show = false)
# plotEarth!(true)

# plot(t,norm.(eachrow(r)).-re)

orbs = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
u = vcat(Vector.(CoeStructure.(orbs)))

plot(t,rad2deg.(scrape.(u,3).-u[1][3]),title="ΔInc vs t (Cowells)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,4).-u[1][4]),title="ΔΩ vs t (Cowells)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,5).-u[1][5]),title="Δω vs t (Cowells)",ylabel="Deg",xlabel="Days",show=true)

RaRp = calcRaRp.(orbs)
plot(t,scrape.(RaRp,1).-re,label="za",show=false)
plot!(t,scrape.(RaRp,2).-re,label="zp",title="za and zp vs time (Cowells)",ylabel="km",xlabel="days",show=true)

println("Deorbit at: $(tf) days")
# RaRp = calcRaRp.(orbs)
# plot(t,scrape.(RaRp,1).-re,label="ra")
# plot!(t,scrape.(RaRp,2).-re,label="rp",title="ra and rp vs time")

# rai,ra = findmaxima(norm.(eachrow(r)).-re)
# rpi,rp = findminima(norm.(eachrow(r)).-re)
# rat = t[rai]
# rpt = t[rpi]
# plot(rat,ra)
# plot!(rpt,rp,show=true)


# VOP
println("-VOP-")
# orb₀ = rv_to_kepler(r₀,v₀)

X₀ = Vector(CoeStructure(orb₀))

function EOM(X,p,t)
    X = CoeStructure(X...) # Convert to struct
    μ = p.μ
    BC = p.BC
    r,v = kepler_to_rv(KeplerianElements(X))
    r/=1000
    v/=1000
    Ẋ =  Vector(calcGaussVOP(X,calcDrag(Vector(r),Vector(v),BC))) # convert to Vector
end

function myCondition(X,t,integrator)
    X = CoeStructure(X...) # Convert to struct
    (X.a-re)<100
end

myCB = DiscreteCallback(myCondition,terminate!)

tspan = (0,1000*24*3600)
myProb = ODEProblem(EOM,X₀,tspan,p)
sol = @time solve(myProb,reltol = 1e-12,abstol = 1e-12,callback = myCB)

t = sol.t
t /= (24*3600)
tf = t[end]

u = sol.u
orbs = KeplerianElements.(CoeStructure.(u))

plot(t,rad2deg.(scrape.(u,3).-u[1][3]),title="ΔInc vs t (VOP)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,4).-u[1][4]),title="ΔΩ vs t (VOP)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,5).-u[1][5]),title="Δω vs t (VOP)",ylabel="Deg",xlabel="Days",show=true)

RaRp = calcRaRp.(orbs)
plot(t,scrape.(RaRp,1).-re,label="za",show=false)
plot!(t,scrape.(RaRp,2).-re,label="zp",title="za and zp vs time (VOP)",ylabel="km",xlabel="days",show=true)

println("Deorbit at: $(tf) days")
# Enkes
println("-Enkes-")
Δt = 100#s
X₀ = [r₀;v₀] |> Vector
X = [X₀]
t = 0
BC = p.BC
δv = [0,0,0]
δr = [0,0,0]

function calcδa(ap,δr,r,rosc)
q = dot(δr,(r.*2-δr))/norm(r)^2
F = q*(q^2-3q+3)/(1+(1-q)^(3/2))
ap-(μe/rosc^3)*(δr-F*r)
end

@time while (norm(r)-re)>100||t>120*24*3600
    rectify = true
    state = X[end]
    global r = state[1:3]
    global v = state[4:6]
    sv = OrbitStateVector(t,r*1000,v*1000)
    prop = Propagators.init(Val(:TwoBody),sv_to_kepler(sv))
    rosc,vosc = Propagators.step!(prop, Δt)./1000
    ap = calcDrag(r,v,BC)
    if rectify
        δa = ap
    else
        δa=calcδa(ap,δr,r,rosc)
    end
    # δa = [0,0,0]
    global δv = δa*Δt+δv
    global δr = .5δa*Δt^2+δv*Δt+δr
    if rectify
        r = δr+rosc
        v = δv+vosc
        δr = [0,0,0]
        δv = [0,0,0]
    end
    state = vcat(r,v)
    push!(X,state)
    global t += Δt
    # println("loop t: $(t) rectify: $(rectify)")
end

r = hcat(scrape.(X,1),scrape.(X,2),scrape.(X,3))
v = hcat(scrape.(X,4),scrape.(X,5),scrape.(X,6))
t = 0:Δt:t
t /= (24*3600)
tf = t[end]
orbs = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
u = vcat(Vector.(CoeStructure.(orbs)))

plot(t,rad2deg.(scrape.(u,3).-u[1][3]),title="ΔInc vs t (Enkes)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,4).-u[1][4]),title="ΔΩ vs t (Enkes)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,5).-u[1][5]),title="Δω vs t (Enkes)",ylabel="Deg",xlabel="Days",show=true)

RaRp = calcRaRp.(orbs)
plot(t,scrape.(RaRp,1).-re,label="za",show=false)
plot!(t,scrape.(RaRp,2).-re,label="zp",title="za and zp vs time (Enkes)",ylabel="km",xlabel="days",show=true)

println("Deorbit at: $(tf) days")

# 2 body
println("-2 body-")
t = 0:100:120*24*3600

prop = Propagators.init(Val(:TwoBody),orb₀)
ret = @time Propagators.propagate!.(prop, t)
r = first.(ret)
v = last.(ret)
v./=1000
r./=1000

t = t|>Vector
t /= (24*3600)
tf = t[end]
orbs = Vector{KeplerianElements}

r=transpose(reduce(hcat, r))
v=transpose(reduce(hcat, v))

orbs = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
u = vcat(Vector.(CoeStructure.(orbs)))

plot(t,rad2deg.(scrape.(u,3).-u[1][3]),title="ΔInc vs t (2B)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,4).-u[1][4]),title="ΔΩ vs t (2B)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,5).-u[1][5]),title="Δω vs t (2B)",ylabel="Deg",xlabel="Days",show=true)

RaRp = calcRaRp.(orbs)
plot(t,scrape.(RaRp,1).-re,label="za",show=false)
plot!(t,scrape.(RaRp,2).-re,label="zp",title="za and zp vs time (2B)",ylabel="km",xlabel="days",show=true)

# Problem2
println("--Probelm 2--")
rp = 300 + re
ra = 3092 + re
ecc = (ra-rp)/(ra+rp)
a = (ra+rp)/2

inc = 28 |> deg2rad
Ω = 45|> deg2rad
ω = 30|> deg2rad
θ = 40 |> deg2rad

t = 0

orb₀ = KeplerianElements(t,a*1000,ecc,inc,Ω,ω,θ)
r₀,v₀ = kepler_to_rv(orb₀)
r₀ /=1000 
v₀ /=1000 

X₀ = [r₀;v₀] |> Vector

function EOM1(X,p,t)
    μ = p.μ
    r = X[1:3]
    v = X[4:6]
    Jaccel = calcJaccel(r)
    a = calcAccel(r,μ) + Jaccel[1]
    Ẋ = [v;a]
end
function EOM2(X,p,t)
    μ = p.μ
    r = X[1:3]
    v = X[4:6]
    Jaccel = calcJaccel(r)
    a = calcAccel(r,μ) + Jaccel[1]+Jaccel[2]
    Ẋ = [v;a]
end

tspan = (0,48*3600)
myProb = ODEProblem(EOM1,X₀,tspan,p)
sol = solve(myProb,reltol = 1e-12,abstol = 1e-12)
r = [sol[1,:] sol[2,:] sol[3,:]]
v = [sol[4,:] sol[5,:] sol[6,:]]
t1 = sol.t
t1 /= (24*3600)
# tf = t[end]
# plot(eachcol(r)...,title = "48hrs J2",show = false)
# plotEarth!(true)
orbs1 = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
u1 = vcat(Vector.(CoeStructure.(orbs1)))

myProb = ODEProblem(EOM2,X₀,tspan,p)
sol = solve(myProb,reltol = 1e-12,abstol = 1e-12)
r = [sol[1,:] sol[2,:] sol[3,:]]
v = [sol[4,:] sol[5,:] sol[6,:]]
t2 = sol.t
t2 /= (24*3600)
# tf = t[end]
# plot(eachcol(r)...,title = "48hrs J2",show = false)
# plotEarth!(true)
orbs2 = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
u2 = vcat(Vector.(CoeStructure.(orbs2)))

plot(t1,rad2deg.(scrape.(u1,3).-u1[1][3]),title="ΔInc vs t (48hrs J2)",ylabel="Deg",xlabel="Days",show=false,label = "J2")
plot!(t2,rad2deg.(scrape.(u2,3).-u2[1][3]),title="ΔInc vs t (48hrs J2)",ylabel="Deg",xlabel="Days",show=true,label = "J3")
plot(t1,rad2deg.(scrape.(u1,4).-u1[1][4]),title="ΔΩ vs t (48hrs J2)",ylabel="Deg",xlabel="Days",show=false,label = "J2")
plot!(t2,rad2deg.(scrape.(u2,4).-u2[1][4]),title="ΔΩ vs t (48hrs J2)",ylabel="Deg",xlabel="Days",show=true,label = "J3")
plot(t1,rad2deg.(scrape.(u1,5).-u1[1][5]),title="Δω vs t (48hrs J2)",ylabel="Deg",xlabel="Days",show=false,label = "J2")
plot!(t2,rad2deg.(scrape.(u2,5).-u2[1][5]),title="Δω vs t (48hrs J2)",ylabel="Deg",xlabel="Days",show=true,label = "J3")

RaRp1 = calcRaRp.(orbs1)
RaRp2 = calcRaRp.(orbs2)
plot(t1,scrape.(RaRp1,1).-re,label="za J2",show=false)
plot!(t1,scrape.(RaRp1,2).-re,label="zp J2",title="za and zp vs time (48hrs J2)",ylabel="km",xlabel="days",show=false)
plot!(t2,scrape.(RaRp2,1).-re,label="za J3",show=false)
plot!(t2,scrape.(RaRp2,2).-re,label="zp J3",title="za and zp vs time (48hrs J2)",ylabel="km",xlabel="days",show=true)


#Problem 3
println("--Probelm 3--")
SpaceIndices.init()
jd = DateTime(Date(2000,11,18),Time(4,20,0))
p = Params(p.μ,p.BC,jd)

function EOM(X,p,t)
    μ = p.μ
    BC = p.BC
    r = X[1:3]
    v = X[4:6]
    a = calcAccel(r,μ) + calcDrag(r,v,BC)
    Ẋ = [v;a]
end

function myCondition(X,t,integrator)
    X = CoeStructure(X...) # Convert to struct
    (X.a-re)<100
end

myCB = DiscreteCallback(myCondition,terminate!)

m = 100
Cd = 2.2
A = π*(1/2)^2
BC = m/(Cd*A)

rp = 215 + re
ra = 939 + re
ecc = (ra-rp)/(ra+rp)
a = (ra+rp)/2

inc = 65.2 |> deg2rad
Ω = 340|> deg2rad
ω = 58|> deg2rad
θ = 332 |> deg2rad

t = 0

orb₀ = KeplerianElements(t,a*1000,ecc,inc,Ω,ω,θ)
r₀,v₀ = kepler_to_rv(orb₀)
r₀ /=1000 
v₀ /=1000 

X₀ = [r₀;v₀] |> Vector
function EOM1(X,p,t)
    μ = p.μ
    BC = p.BC
    r = X[1:3]
    v = X[4:6]
    a = calcAccel(r,μ) + calcDrag(r,v,BC)
    Ẋ = [v;a]
end
function EOM2(X,p,t)
    μ = p.μ
    BC = p.BC
    jd = p.jd
    r = X[1:3]
    v = X[4:6]
    h = norm(r)-re
    q=r_eci_to_ecef(Quaternion,J2000(),PEF(),date_to_jd(jd))
    r2 = vecRotate(r,q)
    lat,lon = ecef_to_geodetic(r2)
    atm = AtmosphericModels.nrlmsise00(jd,h*1000,lat,lon)
    ρ = atm.total_density
    a = calcAccel(r,μ) + calcDrag(r,v,BC,ρ=ρ)
    Ẋ = [v;a]
end
function myCondition(X,t,integrator)
    r = X[1:3]
    # v = X[4:6]
    (norm(r)-re)<100
end

myCB = DiscreteCallback(myCondition,terminate!)

println("-Exponential Model-")
tspan = (0,1000*24*3600)
myProb = ODEProblem(EOM1,X₀,tspan,p)
sol = @time solve(myProb,reltol = 1e-8,abstol = 1e-8,callback = myCB)
r = [sol[1,:] sol[2,:] sol[3,:]]
v = [sol[4,:] sol[5,:] sol[6,:]]
t1 = sol.t
t1 /= (24*3600)
tf1 = t1[end]
plot(eachcol(r)...,title = "nrlmsise00",show = false)
plotEarth!(true)

orbs1 = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
u1 = vcat(Vector.(CoeStructure.(orbs1)))

println("-NRLMSISE-")
myProb = ODEProblem(EOM2,X₀,tspan,p)
sol = @time solve(myProb,reltol = 1e-8,abstol = 1e-8,callback = myCB)
r = [sol[1,:] sol[2,:] sol[3,:]]
v = [sol[4,:] sol[5,:] sol[6,:]]
t2 = sol.t
t2 /= (24*3600)
tf2 = t2[end]


orbs2 = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
u2 = vcat(Vector.(CoeStructure.(orbs2)))

plot(t1,rad2deg.(scrape.(u1,3).-u1[1][3]),title="ΔInc vs t (Drag Comparison)",ylabel="Deg",xlabel="Days",show=false,label = "exp")
plot!(t2,rad2deg.(scrape.(u2,3).-u2[1][3]),title="ΔInc vs t (Drag Comparison)",ylabel="Deg",xlabel="Days",show=true,label = "NRL")
plot(t1,rad2deg.(scrape.(u1,4).-u1[1][4]),title="ΔΩ vs t (Drag Comparison)",ylabel="Deg",xlabel="Days",show=false,label = "exp")
plot!(t2,rad2deg.(scrape.(u2,4).-u2[1][4]),title="ΔΩ vs t (Drag Comparison)",ylabel="Deg",xlabel="Days",show=true,label = "NRL")
plot(t1,rad2deg.(scrape.(u1,5).-u1[1][5]),title="Δω vs t (Drag Comparison)",ylabel="Deg",xlabel="Days",show=false,label = "exp")
plot!(t2,rad2deg.(scrape.(u2,5).-u2[1][5]),title="Δω vs t (Drag Comparison)",ylabel="Deg",xlabel="Days",show=true,label = "NRL")

RaRp1 = calcRaRp.(orbs1)
RaRp2 = calcRaRp.(orbs2)
plot(t1,scrape.(RaRp1,1).-re,label="za exp",show=false)
plot!(t1,scrape.(RaRp1,2).-re,label="zp exp",title="za and zp vs time (Drag Comparison)",ylabel="km",xlabel="days",show=false)
plot!(t2,scrape.(RaRp2,1).-re,label="za NRL",show=false)
plot!(t2,scrape.(RaRp2,2).-re,label="zp NRL",title="za and zp vs time (Drag Comparison)",ylabel="km",xlabel="days",show=true)

println("Deorbit time:")
println("Exponential Model: $(tf1) days")
println("NRLMSISE Model: $(tf2) days")
println("End Program")
