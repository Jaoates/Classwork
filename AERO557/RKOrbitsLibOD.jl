module OrbitsLibOD

include("OrbitsLibBase.jl")
using .OrbitsLibBase
using LinearAlgebra, Roots
using SatelliteToolbox

export calcRADec, calcElAz
export lamberts, lamberts_gauss, lamberts_min, gaussIOD, gibbs, herrick_gibbs
export genObsMat!
export Site

struct Site 
    lat::Number #deg
    lon::Number #deg
    h::Number # alt in m
end

function calc_r_ECEF(s::Site)

end


function calcRADec(r_ECI::AbstractVector)
    u_r = r_ECI ./ norm(r_ECI)
    δ = asind(u_r[3])
    if u_r[2] > 0
        α = acosd(u_r[1]/cosd(δ))
    else
        α = 360 - acosd(u_r[1]/cosd(δ))
    end
    return α, δ
end

function calcElAz(r_NED::AbstractVector)
    El = atand(-r_NED[3], sqrt(r_NED[1]^2 + r_NED[2]^2))
    Az = mod(atand(r_NED[2], r_NED[1]), 360)
    return El, Az
end

S(z) = sum(k -> ((-1)^k * z.^k./factorial(2*k+3)), 0:8)
C(z) = sum(k -> ((-1)^k * z.^k./factorial(2*k+2)), 0:8) 

"""
**Gauss initial orbit determination, needs 3 observations** 
*r_site* observation locations in ECI [km]
*ρ_h* unit vec to target in ECI  
*jd* time of observations [julian day]
"""
function gaussIOD(r_site::AbstractVector{T1}, ρ_h::AbstractVector{T2}, jd::AbstractVector; μ=μₑ, ext=true, ext_tol=1e-8) where {T1<:AbstractVector, T2<:AbstractVector}
    length(r_site) == 3 || throw(ArgumentError("Invalid length, should be 3"))
    length(ρ_h) == 3 || throw(ArgumentError("Invalid length, should be 3"))
    length(jd) == 3 || throw(ArgumentError("Invalid length, should be 3"))

    τ1 = (jd[1] - jd[2])*24*3600 # to s
    τ3 = (jd[3] - jd[2])*24*3600

    a1 = τ3/(τ3 - τ1)
    a1u = τ3*((τ3 - τ1)^2 - τ3^2) / (6*(τ3 - τ1))
    a3 = -τ1 / (τ3 - τ1)
    a3u = -(τ1*((τ3 - τ1)^2 - τ1^2)) / (6*(τ3 - τ1))
    M = inv(reduce(hcat, ρ_h)) * reduce(hcat, r_site) 
    d1 = M[2,1]*a1 - M[2,2] + M[2,3]*a3
    d2 = M[2,1]*a1u + M[2,3]*a3u
    C = dot(ρ_h[2], r_site[2])

    poly(r2) = r2^8 - (d1^2 + 2*C*d1 + norm(r_site[2])^2)*r2^6 - 2*μ*(C*d2 + d1*d2)*r2^3 - μₑ^2*d2^2 
    r2n = find_zero(poly, (0, 1e15))
    u = μ/r2n^3
    c1 = a1 + a1u*u
    c3 = a3 + a3u*u
    c2 = -1

    calcρ(c::AbstractVector, M::AbstractMatrix) = (M * -c) ./ c
    ρ = calcρ([c1, c2, c3], M)
    r = ρ.*ρ_h + r_site
    
    f1 = 1 - τ1^2 * .5*μ/r2n^3
    f3 = 1 - τ3^2 * .5*μ/r2n^3
    g1 =  τ1 - τ1^3*(1/6)*μ/r2n^3
    g3 =  τ3 - τ3^3*(1/6)*μ/r2n^3
    v2 = (-f3*r[1] + f1*r[3]) / (f1*g3 - f3*g1)
    if ext
        ρ_old = zeros(3)
        while norm(ρ - ρ_old) > ext_tol 
            Δθ_1 = acosd(dot(r[1], r[2])/(norm(r[1]) * norm(r[2])))
            Δθ_3 = acosd(dot(r[3], r[2])/(norm(r[3]) * norm(r[2])))
            h = cross(r[2],v2)
            p = norm(h)^2/μ
            f1 = 1 - (norm(r[1])/p)*(1 - cosd(Δθ_1))
            f3 = 1 - (norm(r[3])/p)*(1 - cosd(Δθ_3))
            g1 = -norm(r[1])*norm(r[2])*sind(Δθ_1)/sqrt(μ*p)
            g3 = norm(r[3])*norm(r[2])*sind(Δθ_3)/sqrt(μ*p)
            c1 = g3/(f1*g3 - f3*g1)
            c3 = -g1/(f1*g3 - f3*g1)
            ρ_old = ρ
            ρ = calcρ([c1, c2, c3], M)
            r = ρ.*ρ_h + r_site
            v2 = (-f3*r[1] + f1*r[3]) / (f1*g3 - f3*g1)
        end
    end
    return v2, r[2]
end

function lamberts(r1::AbstractVector, r2::AbstractVector, Δt::Float64; tm::Int=1, tol = 1e-8, μ = μₑ)
    @assert tm == 1 || tm == -1
    cosΔθ = dot(r1, r2)/(norm(r1) * norm(r2))
    Δθ = asind(tm*sqrt(1 - cosΔθ^2))
    
    A = sind(Δθ)*sqrt(norm(r1)*norm(r2) / (1 - cosd(Δθ)))
    @assert abs(A) > tol

    y(z) = norm(r1) + norm(r2) + A*(z*S(z) - 1)/sqrt(C(z))
    F(z) = (y(z)/C(z))^1.5 * S(z) + A * sqrt(y(z)) - sqrt(μ)*Δt
    z_bounds = (-4*pi^2, 4*pi^2)
    z = find_zero(F, 0, rtol=tol)

    f = 1 - y(z)/norm(r1)
    g = A * sqrt(y(z)/μ)
    ḟ = sqrt(μ)/(norm(r1)*norm(r2)) * sqrt(y(z)/C(z)) * (z*S(z) - 1)
    ġ = 1 - y(z)/norm(r2)

    v₁ = 1/g * (r2 - f*r1)
    v₂ = 1/g * (ġ*r2 - r1)
    return v₁, v₂
end

function lamberts_gauss(r1::AbstractVector, r2::AbstractVector, Δt::Float64; tm::Int=1, tol = 1e-6, μ = μₑ)
    @assert tm == 1 || tm == -1
    cosΔθ = dot(r1, r2)/(norm(r1) * norm(r2))
    Δθ = asind(tm*sqrt(1 - cosΔθ^2))

    l = ((norm(r1) + norm(r2)) / (4*sqrt(norm(r1)*norm(r2))*cosd(Δθ/2))) - .5
    m = μ*Δt^2 / (2*sqrt(norm(r1)*norm(r2))*cosd(Δθ/2))^3
    y = 1
    y_old = 0.0
    while (abs(y - y_old) > 1e-6) 
        global x = m/y^2 - l
        X = (4/3)*(1 + (6*x)/5 + (6*8*x^2)/(5*7) + (6*8*10*x^3)/(5*7*9))
        y_old = y
        y = 1 + X*(l+x)
    end
    ΔE = 4*asin(sqrt(x))
    a = ((sqrt(μ)*Δt) / (2*y*sqrt(norm(r1)*norm(r2))*sin(ΔE/2)*cosd(Δθ/2)))^2
    f = 1 - a/norm(r1) * (1 - cos(ΔE))
    g = Δt - a^(1.5)/sqrt(μ) * (ΔE - sin(ΔE))
    v1 = (r2 - f*r1)/g 
    return v1
end

"""
Returns velocity at min energy (semimajor axis) point, times at min energy, and abs min time
"""
function lamberts_min(r1::AbstractVector, r2::AbstractVector; tm::Int=1, μ = μₑ)
    @assert tm == 1 || tm == -1
    cosΔθ = dot(r1, r2)/(norm(r1) * norm(r2))
    Δθ = asind(tm*sqrt(1 - cosΔθ^2))

    c = sqrt(norm(r1)^2 + norm(r2)^2 - 2*norm(r1)*norm(r2)*cosd(Δθ))
    s = (norm(r1) + norm(r2) + c)/2
    a_min = s/2
    p_min = (1 - cosd(Δθ))*norm(r1)*norm(r2)/c
    e_min = sqrt(1 - 2*p_min/s)
    α_e = π
    β_e = 2*asin(sqrt((s - c)/s))
    t_min_e = sqrt(a_min^3/μ)*(α_e - tm*(sin(β_e)))
    t_min_abs = (1/3)*sqrt(2/μ) * (s^(1.5) - (s - c)^(1.5))
    v_1 = (r2 .- (1 - norm(r2)/p_min*(1 - cosd(Δθ)))*r1 ) * sqrt(μ*p_min)/(norm(r1)*norm(r2)*sind(Δθ))
    return v_1, t_min_e, t_min_abs
end

function gibbs(r1::AbstractVector, r2::AbstractVector, r3::AbstractVector;  μ = μₑ)
    α12 = acosd(dot(r1, r2)/(norm(r1) * norm(r2)))
    α23 = acosd(dot(r2, r3)/(norm(r2) * norm(r3)))
     
    if (α12 < 3 || α23 < 3)
        @warn "Separation of $(round(α12, digits=2)) & $(round(α23, digits=2)) degrees, Gibbs may have reduced accuracy"
    end

    z31 = cross(r3, r1)
    z23 = cross(r2, r3)
    if asind(dot(z23, r1)/(norm(z23)*norm(r1))) > 1.0
        @warn "Coplanar-ness above 1 degree, Gibbs may have reduced accuracy"
    end
    z12 = cross(r1, r2)
    D = z12 + z23 + z31
    N = norm(r1)*z23 + norm(r2)*z31 + norm(r3)*z12
    S = r1*(norm(r2)-norm(r3)) + r2*(norm(r3)-norm(r1)) + r3*(norm(r1)-norm(r2))
    v2 = sqrt(μ/(norm(N)*norm(D))) * (cross(D,r2)/norm(r2) + S)
    return v2
end

function herrick_gibbs(r1::AbstractVector, r2::AbstractVector, r3::AbstractVector, JD::AbstractVector;  μ = μₑ)
    α12 = acosd(dot(r1, r2)/(norm(r1) * norm(r2)))
    α23 = acosd(dot(r2, r3)/(norm(r2) * norm(r3)))
     
    if (α12 > 1 || α23 > 1)
        @warn "Separation of $(round(α12, digits=2)) & $(round(α23, digits=2)) degrees, Herrick Gibbs may have reduced accuracy"
    end

    Δt31 = (JD[3] - JD[1])*24*3600
    Δt21 = (JD[2] - JD[1])*24*3600
    Δt32 = (JD[3] - JD[2])*24*3600
    v2 = -Δt32 * ( (Δt31*Δt21)^-1 + μ/(12*norm(r1)^3)) * r1 +
        (Δt32 - Δt21) * ((Δt32*Δt21)^-1 + μ/(12*norm(r2)^3)) * r2 +
        Δt21 * ((Δt32*Δt31)^-1 + μ/(12*norm(r3)^3)) * r3

    return v2
end
"""
Returns m x n matrix of maximally spaced elements of obsVec, where m is length of obsVec ÷ n
If the length of obsVec is not divisible by n, middle elements will be discarded until it is, modifying the input
"""
function genObsMat!(obsVec::AbstractVector; n::Integer = 3)
    N = length(obsVec)
    while N % n != 0
        deleteat!(obsVec, N÷2)
        N = length(obsVec)
    end
    m::Integer = N / 3
    return reduce(hcat, [obsVec[1+(i*m):m+(i*m)] for i in 0:n-1])
end

end