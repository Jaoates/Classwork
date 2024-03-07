# AERO 452 HW1
# Problem 1
include("JoshOrbitsLib.jl")
using SatelliteToolbox
using Plots
using ReferenceFrameRotations
using DifferentialEquations

plotlyjs()

# show?
showTrue = true

# assume an initial epoch
T0 = date_to_jd(-4713, 11, 24, 12, 0, 0.0)

h = 51400; # km2/s
e = 0.0006387;
i = 51.65 |> deg2rad; #degs
Ω = 15 |> deg2rad; #degs
ω = 157 |> deg2rad; #degs
θ = 15 |> deg2rad; #degs
a = calcSemimajor(h,e)
a *= 1000

A = KeplerianElements(T0,a,e,i,Ω,ω,θ)
Pa = orbital_period(A,perturbation=:J0)
Ap = Propagators.init(Val(:TwoBody), A)

numPeriods = 10
Tspan = 0:1:Pa*numPeriods
# Tspan = 0:1:10000

h = 51398; # km2/s
e = 0.0072696;
i = 50 |> deg2rad; #degs
Ω = 15 |> deg2rad; #degs
ω = 140 |> deg2rad; #degs
θ = 15 |> deg2rad; #degs
a = calcSemimajor(h,e)
a *= 1000

B = KeplerianElements(T0,a,e,i,Ω,ω,θ)
Bp = Propagators.init(Val(:TwoBody), B)

ret = Propagators.propagate!.(Ap, Tspan)
ra = first.(ret)
va = last.(ret)
va = va./1000
ra = ra./1000

ret = Propagators.propagate!.(Bp, Tspan)
rb = first.(ret)
vb = last.(ret)
vb = vb./1000
rb = rb./1000

ra=transpose(reduce(hcat, ra))
rb=transpose(reduce(hcat, rb))


va=transpose(reduce(hcat, va))
vb=transpose(reduce(hcat, vb))

# plot(ra[:,1],ra[:,2],ra[:,3])
# plot!(rb[:,1],rb[:,2],rb[:,3])
plot(eachcol(ra)...,show = false)
plot!(eachcol(rb)...,show = false)

dr = rb-ra
# plot(dr[:,1],dr[:,2],dr[:,3],show = showTrue)
plot(eachcol(dr)...,show = false)

B_ECI = eye(3)


rRel = zeros(length(Tspan),3)
for (i,t) in enumerate(Tspan)
    rai = ra[i,:]
    rbi = rb[i,:]
    vai = va[i,:]
    vbi = vb[i,:]
    dri = dr[i,:]

    B_LVLH = calcLVLHBasis(rai,vai)
    Q_ECI2LVLH = vectrixDot(B_ECI,B_LVLH)
    rRel[i,:] = vecRotate(dr[i,:],Q_ECI2LVLH)
end

dist = norm.(eachrow(rRel))

plot(eachcol(rRel)...,show = false)
plot(Tspan/3600,dist,show = showTrue,title = "Distance Over Time of s/c",xaxis = "Time [hr]",yaxis = "Distance [km]")

(m,i) = findmin(dist)
scatter!([Tspan[i]/3600],[m],label="Closest Approach")


println("PROBLEM 1")
println("The closest approach in km is at a distance of:")
println(m)
println("and in hours since T0 this occurs:")
println(Tspan[i]/3600)

#Problem 2

# Spacecraft A and B are in coplanar, circular geocentric orbits. The orbital radii are shown in the
# figure. When B is directly below A, as shown, calculate B’s acceleration arel)xyz relative to A.
za = 300
zb = 250
re = EARTH_EQUATORIAL_RADIUS/1000
ra = (za+re)*[1,0,0]
rb = (zb+re)*[0,1,0]

# ra = 8000 * [1,0,0]  #km
# rb = 7000 * [1,0,0]

va = calcSpeedCirc(norm(ra)) * [0,1,0]
vb = calcSpeedCirc(norm(rb)) * [-1,0,0]

(rRel,vRel,aRel) = calcRelativeKinematics(ra,va,rb,vb)

B_ECI = eye(3)
B_LVLH = calcLVLHBasis(ra,va)
Q_ECI2LVLH = vectrixDot(B_ECI,B_LVLH)
rRel = vecRotate(rRel,Q_ECI2LVLH)
vRel = vecRotate(vRel,Q_ECI2LVLH)
aRel = vecRotate(aRel,Q_ECI2LVLH)

println("PROBLEM 2")
println("In order and in km, km/s, and km/s^2 respectively, are the relative position, velocity, and acceleration of the s/c in LVLH below:")
println(rRel)
println(vRel)
println(aRel)

# Problem 3

zp = 250
rp = zp + EARTH_EQUATORIAL_RADIUS/1000

T0 = date_to_jd(-4713, 11, 24, 12, 0, 0.0)

e = 0.1;
i = 51 |> deg2rad; #degs
Ω = 0 |> deg2rad; #degs
ω = 0 |> deg2rad; #degs
θ = 0 |> deg2rad; #degs

a = calcSemimajor(rp,e,input = "rpe")*1000

A = KeplerianElements(T0,a,e,i,Ω,ω,θ)
T = calcPeriod(A)
A = kepler_to_rv(A)
T = 90*60
r0 = A[1]/1000
v0 = A[2]/1000
h = norm(cross(r0,v0))

delr0 = [-1, -1, 0] #km/s 
delv0 = 1e-3*[0, 2, 0] # km/s

# function linearRelatieEOM(dX, X, p, t)
#     Dx,Dy,Dz,dDx,dDy,dDz = X[1:6]
#     μ,h = p
#     R̂ = 
#     R = norm(R̂)
#     ddDx = ((2*μ/R^3)+(h^2/R^2))*Dx - ((2*dot(V,R))*h/R^4)*Dy + (2*h/R^2)*dDy
#     ddDx = ((h^2/R^2)-(μ/R^3))*Dy - ((2*dot(V,R))*h/R^4)*Dx - (2*h/R^2)*dDx
#     ddDz = -(μ/R^3)*Dz
#     δ̇r = [dDx,dDy,dDz]
#     δ̈r = [ddDx,ddDy,ddDz]
# end

function linearRelativeEOM(dX, X, p, t)
    # Dx,Dy,Dz,dDx,dDy,dDz = X[1:6]
    # R̂ = X[7:9]
    # V̂ = X[9:12]

    μ, h = p
    dDx,dDy,dDz = dX[:, 1]
    V̂ = dX[:, 2]
    Dx,Dy,Dz = X[:, 1]
    R̂ = X[:, 2]
    # print(V̂)
    # println(R̂)

    R = norm(R̂)
    ddDx = ((2*μ/R^3)+(h^2/R^4))*Dx - ((2*dot(V̂,R̂))*h/R^4)*Dy + (2*h/R^2)*dDy
    ddDy = ((h^2/R^4)-(μ/R^3))*Dy + ((2*dot(V̂,R̂))*h/R^4)*Dx - (2*h/R^2)*dDx
    ddDz = -(μ/R^3)*Dz
    δ̇r = [dDx,dDy,dDz]
    δ̈r = [ddDx,ddDy,ddDz]

    R̈̂ = calcAccel(R̂)
    # print(δ̈r)
    # println(R̈̂)
    return[δ̈r R̈̂]

end

n = 10 # number of periods

p = (μe, h)
tspan = (0, n*T)
prob = SecondOrderODEProblem(linearRelativeEOM, [delv0 v0], [delr0 r0], tspan, p)
sol = solve(prob, saveat=.5, reltol=1e-8)
delr = sol[7:9, :]
plot(delr[2,:],delr[1,:],show = showTrue,xaxis = "y [km]",yaxis = "x [km]",title = "Relative Motion for 10 Orbits")

println("PROBLEM 3")
println("The final seperation of the objects in LVLH in km is:")
println(delr[:,end])
println("The norm of which is:")
println(norm(delr[:,end]))


# Problem 4
T = 90
T = T*60
n = 2π/T

tf = (1/4)*T

δr₀ = 1*[1,0,0]
δ̇r₀ = 10*[0,1,0] / 1000
δr,δṙ= cwEquations(δr₀,δ̇r₀,tf,n)

println("PROBLEM 4")
println("The final distance in 15 min in km based on the CW equations is:")
println(norm(δr))
