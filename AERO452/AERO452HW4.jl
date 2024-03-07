include("JoshOrbitsLib.jl")
# include("OtherOrbitsLib.jl")
using SatelliteToolbox
using Plots
# using ReferenceFrameRotations
using DifferentialEquations
# using Revise

# using Symbolics
# # Symbolics toolbox cites the following:
# #     @article{gowda2021high,
# #     title={High-performance symbolic-numerics via multiple dispatch},
# #     author={Gowda, Shashi and Ma, Yingbo and Cheli, Alessandro and Gwozdz, Maja and Shah, Viral B and Edelman, Alan and Rackauckas, Christopher},
# #     journal={arXiv preprint arXiv:2105.03949},
# #     year={2021}

plotlyjs()

# Problem 1
t = (10,23)
jd = DateTime.(1997,8,t)
jd = date_to_jd.(jd)
nPoints =100000
jd = Vector(range(jd[1],jd[2],nPoints))

sun = sun_position_mod.(jd)./1000
Q = r_eci_to_eci.(Quaternion,MOD(),J2000(),jd)
sun = vecRotate.(sun,Q)

t = range(0,(t[end]-t[1])*24*3600,nPoints) |>Vector

r₀ = [-26175.1034,12757.0706,14626.6556]
v₀ = [2.376641,0.139677,2.078097]
orb₀ = sv_to_kepler(OrbitStateVector(jd[1],r₀*1000,v₀*1000))

prop = Propagators.init(Val(:TwoBody),orb₀)
ret = Propagators.propagate!.(prop, t)

r = first.(ret)./1000
v = last.(ret)./1000
t = t|>Vector

r = map(x->Vector(x),r)
v = map(x->Vector(x),v)

# θA = 90 |> deg2rad
# θB = acos.(re./norm.(r))
# θ = acos.(dot.(sun,r)./(norm.(sun).*norm.(r)))
# F = θA.+θB.>=θ
F = calcInShadow.(r,sun)
F = vcat(0,diff(F))
tstart = t[F .== -1]
tend = t[F .== 1]
tspans = tend.-tstart

eclipseFraction = sum(tspans)/(t[end]-t[1])
println("The fraction of the time in eclipse (and therefor the increase in SRP if shadow was ignored) is: %$(eclipseFraction*100).")
println("The satellite enters and exits eclipse at the following times in days:")
println("eclipse number - start time - end time")
for i = 1:length(tstart)
    println("$(i) - $(tstart[i]/(24*3600)) - $(tend[i]/(24*3600))")
end

# Problem 2

# h = 69,084.1
ecc = 0.741 
Ω = 0|> deg2rad
inc = 63.4|> deg2rad
ω = 270|> deg2rad
θ = 0|> deg2rad
a = 26553.4
T0 = 11.9616
jd₀ = 2454283.0 + T0/24

orb₀ = KeplerianElements(0,a*1000,ecc,inc,Ω,ω,θ)
r₀,v₀ = kepler_to_rv(orb₀)
r₀ /=1000 
v₀ /=1000 
X₀ = [r₀;v₀] |> Vector

struct Params
    μn
    jd₀
end

p = Params(μsun,jd₀)

function EOM(X,p,t)
    μsun = p.μn
    jd₀ = p.jd₀
    r = X[1:3]
    v = X[4:6]
    jd = jd₀ + t/(24*3600)

    rsun = sun_position_mod(jd)/1000
    Q = r_eci_to_eci.(Quaternion,MOD(),J2000(),jd)
    rsun = vecRotate(rsun,Q)

    a = calcAccel(r) + calcNBodyAccel(r,rsun,μsun)
    Ẋ = [v;a]
end

tspan = (0,60*24*3600)
myProb = ODEProblem(EOM,X₀,tspan,p)
sol = solve(myProb,reltol = 1e-12,abstol = 1e-12)
r = [sol[1,:] sol[2,:] sol[3,:]]
v = [sol[4,:] sol[5,:] sol[6,:]]
t = sol.t
t /= (24*3600)
tf = t[end]
orbs = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
u = vcat(Vector.(CoeStructure.(orbs)))

plot(eachcol(r)...,show = true)
plotEarth!(true)

plot(t,rad2deg.(scrape.(u,3).-u[1][3]),title="ΔInc vs t (Solar Gravity)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,4).-u[1][4]),title="ΔΩ vs t (Solar Gravity)",ylabel="Deg",xlabel="Days",show=true)
plot(t,rad2deg.(scrape.(u,5).-u[1][5]),title="Δω vs t (Solar Gravity)",ylabel="Deg",xlabel="Days",show=true)

# RaRp = calcRaRp.(orbs)
# plot(t,scrape.(RaRp,1).-re,label="za",show=false)
# plot!(t,scrape.(RaRp,2).-re,label="zp",title="za and zp vs time (Cowells)",ylabel="km",xlabel="days",show=true)

# Problem 3
r₀ = [
    [-1.05,0,0],
    [-1.15,0,0],
    [.1,0,0],
    [.1,0,0],
    [.1,0,0],
    [1.25,0,0],
    [-.5,0,0],
    [-1.10,0,0]
]
v₀ = [
    [0,-.066429,0],
    [0,-0.0086882909,0],
    [-3.35,3,0],
    [-3.37,3,0],
    [-3.4,3,0],
    [0,0,0],
    [0,0,0],
    [0,0,0]
]
t=[2π,29.46,3.6,6,6,2π,2π,2π]


function EOM(X,p,t)
    r = X[1:3]
    v = X[4:6]
    μs = p
    a = calcCR3BPaccel(r,v,μs)
    Ẋ = [v;a]
end


μs = .01215

X₀ = map((x,y)->[x;y],r₀,v₀ )
tspan = map(x->(0,x),t)
prob = ODEProblem.(EOM,X₀,tspan,μs)
sol =  solve.(prob,reltol = 1e-8,abstol = 1e-8)

r = map(x->x[1:3, :],sol)

titles = range('a','h')
# for i = 1:length(r)
#     plot3d(eachrow(r[i])...,title = titles[i],show = false)
#     plotCanonicalEarthMoon!()
#     plot!(xlims = [-1.5,1],ylims = [-1,1],zlims = [-1,1],show = true)
# end
for i = 1:length(r)
    if titles[i] ∈ ['d','e','g']
        plot!(eachrow(r[i])[1:2]...,show = false,label = "Trajectory $(titles[i])")
    else
        plot(eachrow(r[i])[1:2]...,show = false,label = "Trajectory $(titles[i])")
        scatter!([μs],[0],label = "Earth",show = false)
        scatter!(-[1-μs],[0],label = "Moon",show = false)
    end
    plot!(aspect_ratio = 1,show = true,xlabel = "x [CDU]",ylabel = "y [CDU]")  
end

# Problem 4
n = 500
x =  Vector(range(-2,2,n))
r = Array{Vector}(undef,n,n)
r = [[xi, yi, 0] for xi in x, yi in x]

c = calcZeroVeloJacobi.(r)
x = map(x->x[1],r)
y = map(x->x[2],r)

cutoff = 4
cutoff = c.<cutoff

c = c[cutoff]
x = x[cutoff]
y = y[cutoff]

contour(x,y,c,levels=20,fill = true,show = false)
L = calcLagrangePts()
map(i->scatter!([L[i][1]],[L[i][2]],label="L$(i)"),range(1,5))
plot!(show = true)
