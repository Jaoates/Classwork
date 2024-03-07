using SatelliteToolbox
using LinearAlgebra
using Plots
using ReferenceFrameRotations
using Roots
const μe = GM_EARTH/(1000^3)
const re = EARTH_EQUATORIAL_RADIUS/1000
const ωe = [0,0,1]*72.9211e-6 #rad/s in ECI
const moonDistance = 384400 #km
const moonRadius = 1740 # km
const μsEarthMoon = .01215
const μsun = 132.712*10^9 #km3/s2
const solarConstantEarth = 1376 # W/m^2
const speedOfLight = 299792458 # m/s

function calcSemimajor(p1,p2,μ=μe;input = "he")
    # returns a (semimajor axis) for use with SatelliteToolbox
    if input == "he"
        # takes h (angular momentum) and e (eccentricity) of an orbit, optionally μ (gravitational parameter)
        a = (p1^2/μ)/(1-p2^2);
    elseif input == "rpe"
        # takes radius of perigee (rp) and e (eccentricity) of an orbit, optionally μ (gravitational parameter)
        a = p1/(1-p2)
    else
    print("BadInput")
    end
    return a
end

function calcPeriod(orb::KeplerianElements,μ = GM_EARTH)
    # takes an orbit in KeplerianElements form and gives its period in s
    if orb.e >= 1
        print("Warning: Non-eliptical")
    end
    T = (2π/sqrt(μ))*orb.a^(3/2)
    n = 2π/T
    return (T,n)
end
function calcPeriod(a::Number,μ = μe)
    # takes an orbit in KeplerianElements form and gives its period in s
    T = (2π/sqrt(μ))*a^(3/2)
    n = 2π/T
    return (T,n)
end

function calch(r::Vector,v::Vector)
    h = norm(cross(r,v))
end
function calch(a::Number,ecc::Number;μ::Number=μe)
    h = sqrt(a*(1-ecc^2)*μ)
end
function calch(orb::KeplerianElements;μ::Number = μe)
    h = cach(a,ecc,μ)
end

function vecRotate(v,R::DCM)
    R*v
end

function vecRotate(v,R::Matrix)
    R*v
end

function vecRotate(v,R::Quaternion)
    v  = vect(R \ v * R)
end

function vectrixDot(F1::Matrix,F2::Matrix)
    # % takes the vectrix dot product where the rotation matrix C21 =
    # % vectrixDot(F1,F2) such that r represented in F2 can be found from r2
    # % = C21 * r1. F1 and F2 are 3x3 matricies with columns [x,y,z] where x y
    # % and z are column vectors representing the basis vectors of F1 and F2
    
    x1 = F1[:,1];
    y1 = F1[:,2];
    z1 = F1[:,3];
    
    x2 = F2[:,1];
    y2 = F2[:,2];
    z2 = F2[:,3];
    
    C = [
        dot(x2,x1) dot(x2,y1) dot(x2,z1)
        dot(y2,x1) dot(y2,y1) dot(y2,z1)
        dot(z2,x1) dot(z2,y1) dot(z2,z1)
        ];
end

function eye(n)
    Matrix{Float64}(I(n))
end

function calcLVLHBasis(r::Vector,v::Vector)
    # the lvlh basis in ECI
    x = r/norm(r)
    h = cross(r,v)
    z = h/norm(h)
    y = cross(z,x)
    B_LVLH = [x y z]
end

function calcSpeedCirc(r::Number,μ=μe)
    # r is radius of circular orbit in km
    # v is in km/s
    v = sqrt(μ/r)
end

function calcAngularVelo(r::Vector,v::Vector)
    # r is the position vector of the sc in ECI
    # v is the velocity vector of the sc in ECI
    # returns angular velocity(Ω) in ECI and its derivative (dΩ)
    Ω = cross(ra,va)/(norm(ra)^2)
    dΩ = (-2*dot(va,ra)*Ω)/norm(ra)
    return (Ω,dΩ)
end

function calcAccel(r::AbstractVector,μ=μe)
    # calculates the acceleration of an orbit in free fall as a vector
    accel = -r*(μ/(norm(r)^3))
end

function calcAccel(r::Number,μ=μe)
    # calculates the acceleration of an orbit in free fall as a number
    accel = μ/r^2
end

function calcRelativeKinematics(ra::Vector,va::Vector,rb::Vector,vb::Vector)
    rRel = (rb-ra) 

    # Ω = cross(ra,va)/norm(ra)^2
    # dΩ = -2*dot(va,ra)*Ω/norm(ra)
    (Ω,dΩ) = calcAngularVelo(ra,va)
    
    aa = calcAccel(ra)
    ab = calcAccel(rb)
    
    vRel = vb - va - cross(Ω,rRel)
    aRel = ab - aa - cross(dΩ,rRel) - cross(Ω,cross(Ω,rRel))-2*cross(Ω,vRel)
    return (rRel,vRel,aRel)
end

function cwCoefs(t,n)
    Φrr = [
        4-3*cos(n*t)     0 0
        6*(sin(n*t)-n*t) 1 0
        0                0 cos(n*t)
    ]
    Φrv = [
        (1/n)*sin(n*t)      (2/n)*(1-cos(n*t))          0
        (2/n)*(cos(n*t)-1)  (1/n)*(4*sin(n*t)-3*n*t)    0
        0                   0                           (1/n)*sin(n*t)
    ]
    Φvr = [
        3*n*sin(n*t)        0   0
        6*n*(cos(n*t)-1)    0   0
        0                   0   -n*sin(n*t)
    ]
    Φvv = [
        cos(n*t)    2*sin(n*t)      0
        -2*sin(n*t) 4*cos(n*t)-3    0
        0           0               cos(n*t)
    ]
    return (Φrr,Φrv,Φvr,Φvv)
end

function cwEquations(δr₀::Vector,δ̇r₀::Vector,t::Number,n::Number)
    Φrr,Φrv,Φvr,Φvv = cwCoefs(t,n)
    δr = Φrr * δr₀ + Φrv * δ̇r₀
    δṙ = Φvr * δr₀ + Φvv * δ̇r₀
    return (δr,δṙ)
end

function plotEarth!(show = false)
    n = 100
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    x = cos.(u) * sin.(v)'
    y = sin.(u) * sin.(v)'
    z = ones(n) * cos.(v)'
    x*=re
    y*=re
    z*=re
    surface!(x, y, z;colorbar = false,c = :blues,show = show)
end

function calcDrag(
    r::Vector,
    v::Vector,
    BC::Number;
    ρ::Number = NaN,
    ω::Vector = ωe
    )
    ## the typical drag on a spacecraft due to atmosphere
    # r is position of s/c in ECI in [km]
    # v is velocity of s/c in ECI in [km/s]
    # BC is balistic coeff ( m/(Cd*A) ) in [kg/m^2]
    # ρ is density in [kg/m^3]
    
    if(isnan(ρ))
        # ρ = AtmosphericModels.exponential((norm(r)-re)*1000) # SatelliteToolbox
        ρ = atmosphere((norm(r)-re)) # Curtis
    end
    vrel = v-cross(ω,r)
    BC *= (1000^2)
    ρ *= (1000^3)
    p = (-.5*ρ*norm(vrel)*vrel/BC)
end

function calcEccA(ra::Number,rp::Number)
    ecc = (ra-rp)/(ra+rp)
    a = (ra+rp)/2
    return (ecc,a)
end
function calcEccA(ra::Vector,rp::Vector)
    calcEccA(norm(ra),norm(rp))
end

function calcRaRp(ecc::Number,a::Number)
    # ecc = (ra-rp)/(ra+rp)
    # a = (ra+rp)/2
    ra = a*(ecc + 1)
    rp = a - a*ecc
    return [ra,rp]
end
function calcRaRp(orb::KeplerianElements)
    calcRaRp(orb.e,(orb.a)./1000)
end


function scrape(A::AbstractVector,col::Integer = 1)
    # function allows you to broadcast over a nested array and scrape off a given indicie of each element
    # ex shows how to turn a vector of vectors to a matrix 
    # A = hcat(scrape.(A,1),scrape.(A,2))
    A[col]
end
function scrape(A::Tuple,col::Integer = 1)
    A[col]
end
# function Base.Matrix(A::Vector)
#     mapreduce(permutedims, vcat, A)    
# end
function Base.convert(::Type{Matrix},A::Vector)
    mapreduce(permutedims, vcat, A)
end

struct CoeStructure
    # NO guarantee this thing is ACTUALLY a coe. It just has the structure of one.
    # used to hold time derivatives of coes too. It might make no sense to cast it to an orbit type.
    # one example of this is that addition is defined on it.
    a::Number
    ecc::Number
    inc::Number
    Ω::Number
    ω::Number
    θ::Number
end
function CoeStructure(orb::KeplerianElements)
    CoeStructure(orb.a/1000,orb.e,orb.i,orb.Ω,orb.ω,orb.f)
end
function CoeStructure(orb::Vector)
    CoeStructure(orb...)
end
function SatelliteToolbox.KeplerianElements(orb::CoeStructure; t::Number=0)
    println("running")
    SatelliteToolbox.KeplerianElements(t,orb.a*1000,orb.ecc,orb.inc,orb.Ω,orb.ω,orb.θ)
end
function Base.convert(::Type{T},orb::CoeStructure) where {T<:AbstractVector}
    T([orb.a,orb.ecc,orb.inc,orb.Ω,orb.ω,orb.θ])
end
function Base.Vector(orb::CoeStructure)
    convert(Vector,orb)
end
function Base.:+(o1::CoeStructure,o2::CoeStructure)
    CoeStructure(o1.a+o2.a, o1.e+o2.e, o1.i+o2.i, o1.Ω+o2.Ω, o1.ω+o2.ω, o1.θ+o2.θ)    
end
function Base.:*(o1::CoeStructure,x::Number)
    CoeStructure(o1.a*x, o1.e*x, o1.i*x, o1.Ω*x, o1.ω*x, o1.θ*x)    
end
# function calcMeanCoeStructure(Coes::Vector{CoeStructure})
    #     # function which takes the mean of a CoeStructure, Its a bit funky...

    #     n = length(Coes)
    #     a = sum([c.a for c in Coes])/n
    #     ecc = sum([c.ecc for c in Coes])/n
    #     inc = sum([c.inc for c in Coes])/n
    #     Ω = sum([c.Ω for c in Coes])/n
    #     ω = sum([c.ω for c in Coes])/n
    #     θ = sum([c.θ for c in Coes])/n

    #     return CoeStructure(a,ecc,inc,Ω,ω,θ)
# end


function calcGaussVOP(orb::CoeStructure,p::Vector;μ = μe)
    # orb is the current coes
    # p is the perturbational acceleration in ECI

    # rotate p into LVLH instead of ECI
    B_ECI = eye(3)
    r,v = kepler_to_rv(KeplerianElements(orb))
    r/=1000
    v/=1000
    r = Vector(r)
    v = Vector(v)
    B_LVLH = calcLVLHBasis(r,v)
    Q_ECI2LVLH = vectrixDot(B_ECI,B_LVLH)
    p = vecRotate(p,Q_ECI2LVLH)

    a = orb.a
    ecc = orb.ecc
    inc = orb.inc
    Ω = orb.Ω
    ω = orb.ω
    θ = orb.θ

    pr,ps,pw = p
    h = calch(a,ecc)
    # r1 =h^2/(μ*(1 + ecc*cos(θ)))
    r = norm(r)
    # println(r1-r)

    dh = r*ps
    decc = (h/μ) * sin(θ) * pr + (ps/(μ*h))*( (h^2+μ*r)*cos(θ) + μ*ecc*r )
    dθ = h/r^2 + (1/(ecc*h)) * ( (h^2/μ)*cos(θ)*pr - (r+h^2/μ)*sin(θ)*ps )
    dΩ = r/(h*sin(inc))*sin(ω+θ)*pw
    dinc = (r/h)*cos(ω+θ)*pw
    dω = -1/(ecc*h)*((h^2/μ)*cos(θ)*pr - (r+h^2/μ)*sin(θ)*ps) - pw*r*sin(ω+θ)/(h*tan(inc))

    da = (2*h^2*ecc*decc)/(μ*(ecc^2 - 1)^2) - (2*h*dh)/(μ*(ecc^2 - 1))

    CoeStructure(da,decc,dinc,dΩ,dω,dθ)
end

function calcJaccel(r;re = re,μ = μe)
    #calculate the first six terms of the gravitational perturbations of earth
    #static model,  no options for  J parameters
 
    ri,rj,rk = r
    R = re
    rn = norm(r)
    #J2
    J2 = 1.08262668355e-3
    a2 = -r.*(3*J2*μ*R^2/(2*rn^5)).*([1,1,3].-(5*rk^2/rn^2))
    #J3
    J3 = -2.53265648533e-6
    a3 = -[ri,rj,1].*(5*J3*μ*R^3/(2*rn^7)).*([3,3,6*rk].*rk-[1,1,rk].*(7*rk^3/rn^2)-[0,0,(3/5)*rn^2])
    #J4
    J4 = -1.61962159137e-6
    a4 = r.*(15*J4*μ*R^4/(8*rn^7)).*([1,1,5].-[14,14,70/3].*(rk^2/rn^2).+21*rk^4/rn^4)
    #J5
    J5 = -2.27296082869e-7
    a5 = r.*(3*J5*μ*R^5*rk/(8*rn^9)).*([35,35,105].-[210,210,315].*(rk^2/rn^2).+231*rk^4/rn^4)
    #J6
    J6 = 5.40681239107e-7
    a6 = -r.*(J6*μ*R^6/(16*rn^9)).*([35,35,245].-[945,945,2205].*(rk^2/rn^2).+[3465,3465,4851].*(rk^4/rn^4).+3003*rk^6/rn^6)
    [a2,a3,a4,a5,a6]
end

function calcCR3BPaccel(r::Vector,v::Vector,μs::Number=μsEarthMoon)
    #r is position of the s/c in canonical units
    #v is velocity of the s/c in canonical units
    # the reference frame is the non-inertial berrycentric frame
    # x points to the primary body
    # the 3B system has angular velocity in the z direction
    # y completes

    x,y,z=r
    ẋ,ẏ,ż=v
    r1 = sqrt((x-μs)^2+y^2+z^2)
    r2 = sqrt((x+1-μs)^2+y^2+z^2)
    ẍ = x + 2ẏ - (1-μs)*(x-μs)/r1^3 - μs*(x+1-μs)/r2^3
    ÿ = y - 2ẋ - (1-μs)*y/r1^3 - μs*y/r2^3
    z̈ = -(1-μs)*z/r1^3 - μs*z/r2^3
    a = [ẍ,ÿ,z̈]
    return a
end

function plotSphere!(r,radius,show = false,c = :blues)
    n = 100
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    x = cos.(u) * sin.(v)'
    y = sin.(u) * sin.(v)'
    z = ones(n) * cos.(v)'
    x*=radius
    y*=radius
    z*=radius

    surface!(x.+r[1], y.+r[2], z.+r[3];colorbar = false,c,show = show)
end

function plotCanonicalEarthMoon!(show = false)
    plotSphere!([μs,0,0],re/moonDistance,show)
    plotSphere!([1-μs,0,0],moonRadius/moonDistance,show,:bilbao)
end

function calcZeroVeloJacobi(r::Vector,μs::Number = μsEarthMoon)
    # r is position in canonical units
    # the reference frame is the non-inertial berrycentric frame
    # x points to the primary body
    # the 3B system has angular velocity in the z direction
    # y completes

    x,y,z=r
    r1 = sqrt((x-μs)^2+y^2+z^2)
    r2 = sqrt((x+1-μs)^2+y^2+z^2)
    c = (x^2+y^2)+2*(1-μs)/r1 + 2μs/r2
end

function calcLagrangePts(μs::Number = μsEarthMoon)
    rb1 = μs
    rb2 = 1-μs
    f(x) = x - sign(x-μs)*(1-μs)/(x-μs)^2 - sign((x+1-μs))*μs/(x+1-μs)^2
    intervals = [(-rb2,rb1),(rb1,10),(-10,-rb2)]
    # z = [find_zeros(f,-rb2,rb1),find_zeros(f,rb1,10),find_zeros(f,-10,-rb2)]
    L1,L2,L3 = [append!(find_zeros(f,i),[0,0]) for i in intervals]
    L4 = [μs-.5,sqrt(3)/2,0]
    L5 = [μs-.5,-sqrt(3)/2,0]
    return(L1,L2,L3,L4,L5)
    
#
    # myFun1(x) = x - (1-μs)*(x-μs)/(x-μs)^3 - μs*(x+1-μs)/(x+1-μs)^3
    # myFun3(x) = x - (1-μs)/(x-μs)^2 - μs/(x+1-μs)^2
    # myFun4(x) = (1-μs)*(x+μs)/abs(x+μs)^3 + (x+μs-1)/abs(x+μs-1)^3-x 



    # x = -2:.01:2
    # using Plots
    # plotlyjs()
    # plot(x,myFun1.(x),ylim = (-10,10))
    # plot!(x,myFun2.(x),ylim = (-10,10))
    # plot!(x,myFun3.(x),ylim = (-10,10))
    # plot!(x,myFun4.(x),ylim = (-10,10))
end

function calcNBodyAccel(r::AbstractVector,rPrimaryNbody::AbstractVector,μNbody::Number)
    # is position of sc in ECI
    # rPrimaryNbody is the position of the Nbody in ECI
    # μ is μ of the primary body
    # μNbody is μ of the Nbody
    rn = rPrimaryNbody
    rns = rn-r
    μn = μNbody
    q = dot(r,(rn.*2-r))/norm(rn)^2
    F = q*(q^2-3q+3)/(1+(1-q)^(3/2))    
    ap = μn/norm(rns)^3 * (F*rn - r)
end

function calcInShadow(r::AbstractVector,rsun::AbstractVector;radius::Number = re)
    # r is the spacecraft location in ECI
    # rsun is the location of the sun in ECI
    # radius is the radius of the body were orbiting
    # computes a simple calc based on the assumption that the sun is an infinitely far point source to tell if the s/c is in shadow.
    # assumes the shadow is a binary property
    if norm(r)<=radius
        throw(DomainError(r,"r is too small, norm(r)<radius"))
    end
    θA = 90 |> deg2rad
    θB = acos(radius/norm(r))
    θ = acos(dot(rsun,r)/(norm(rsun)*norm(r)))
    F = θA+θB>=θ
end

function calcSRPaccel(r::AbstractVector,rsun::AbstractVector,SC::Number;
    radius::Number = re, solarConstant::Number = solarConstantEarth)
    # r is the spacecraft location in ECI
    # rsun is the suns location in ECI
    # SC is Solar Radiation Pressure coeff ( m/(Cr*A) ) in [kg/m^2]  --- (similar to balistic coeff)
        # Cr is a constant between 1 and 2 (typically) and contains information on how much radiation is absorbed by the s/c
    # radius is the radius of the body we're orbiting
    # solarConstant is the Solar constant at the location of the s/c in W/m^2
    # ap is the perturbing acceleration in ECI
    # suffers from the same limitations as calcInShadow
    F= calcInShadow(r,rsun;radius = radius)
    psr = F*(solarConstant/speedOfLight)/SC
    psr /=1000
    ap = -psr*(rsun/norm(rsun))
end

