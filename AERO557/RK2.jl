using Pkg
Pkg.activate("A557")
using Revise
includet("../orbitslib/OrbitsLibOD.jl")
includet("../orbitslib/OrbitsLibBase.jl")
using .OrbitsLibOD
using .OrbitsLibBase
using LinearAlgebra, Statistics, Polynomials
using Plots
using Dates
using CSV, DataFrames
using SatelliteToolbox
using Unzip
using PrettyTables

## Problem 1

xoi = collect(1:10)
yoi = [1, 2, 2, 3, 4, 4, 5, 7, 9, 9]
P = fit(xoi, yoi, 1)
x̂ = P[:]

ϵ = yoi - P.(xoi)
RMS = sqrt(mean(ϵ.^2))
σ = std(ϵ)

# plot(P, extrema(xoi)...)
# scatter!(xoi, yoi)
# Makie.plot(P, extrema(xoi)...)

test = 1.5σ
tvec = ϵ.^2 .< test^2
yoi = yoi[tvec]
xoi = xoi[tvec]

ϵ = yoi - P.(xoi)
RMS = sqrt(mean(ϵ.^2))
σ = std(ϵ)
##

## Problem 2

data = DataFrame(CSV.File("HW 2/prob2.txt"; header=false))
times = DateTime.(data[:,1], dateformat"yyyy mm dd  HH MM SS.s")
JD = date_to_jd.(times)

Az, El, ρ_n = data[:,2], data[:, 3], data[:, end]
bias = [.0224, .0139, 92.5e-3] #az deg, el deg, km
lat, lon, h = [21.5721, -158.2666, 300.2] # deg, deg, m
site = Site(lat, lon, h)

r_site_ECEF = 1e-3*geodetic_to_ecef(site.lat |> deg2rad, site.lon |> deg2rad, site.h) #m

##

# RA, DEC = unzip(AstroLib.hor2eq.(El, Az, JD, site.lat, site.lon, site.h))
# ρ_ECI = ρ_n .* [ [cosd(ra)*cosd(dec), cosd(dec)*sind(ra), sind(dec)] for (ra, dec) in zip(RA, DEC)]
# r_site_ECI = r_ecef_to_eci.(PEF(), J2000(), JD) .* Ref(r_site_ECEF)
# r_ECI = ρ_ECI + r_site_ECI

# ρ_NED = ρ_n .*  [ [cosd(az)*cosd(el), cosd(el)*sind(az), -sind(el)] for (az, el) in zip(Az, El) ]
# ρ_ECEF = 1e-3*ned_to_ecef.(ρ_NED*1e3, site.lat |> deg2rad, site.lon |> deg2rad, site.h)
# r_ECEF = Ref(r_site_ECEF) .+ ρ_ECEF

# r_ECI = r_ecef_to_eci.(PEF(), J2000(), JD) .* r_ECEF

# r_mat = genObsMat!(r_ECI)
# JD_mat = genObsMat!(JD)
# v2 = [herrick_gibbs(vecs..., jds) for (vecs, jds) in zip(eachrow(r_mat), eachrow(JD_mat))]
# SVs = OrbitStateVector.(JD_mat[:,2], 1e3*r_mat[:,2], 1e3*v2)
# r, v = unzip(Propagators.propagate_to_epoch.(Val(:TwoBody), JD[1], sv_to_kepler.(SVs)))

# X_nom_1 = 1e-3*vcat(mean(r), mean(v))

##

struct BatchData
    X_nom::Vector{Float64}
    P::Matrix{Float64}
    HtWH::Matrix{Float64}
    HtWy::Vector{Float64}
    W::Matrix{Float64}
    H::Matrix{Float64}
end

function propagate_state(X_1::AbstractVector, JD::AbstractVector)
    SV_1 = OrbitStateVector(JD[1], 1e3*X_1[1:3], 1e3*X_1[4:6])
    r, v = unzip(Propagators.propagate_to_epoch.(Val(:TwoBody), JD, sv_to_kepler(SV_1)))
    X_i = 1e-3*Matrix(hcat([r'...;], [v'...;]))
    return X_i::AbstractMatrix
end

function state_to_obs(X::AbstractVector, jd::Number, site::Site)
    @assert length(X) == 6 "State must be [rx, ry, rz, vx, vy, vz]"
    r = 1e3*X[1:3] #m
    r_ECEF = r_eci_to_ecef(J2000(), PEF(), jd) * r #m
    r_NED = ecef_to_ned(r_ECEF, site.lat |> deg2rad, site.lon |> deg2rad, site.h; translate=true)
    el, az = calcElAz(r_NED) 
    r_site_ECEF = 1e-3*geodetic_to_ecef(site.lat |> deg2rad, site.lon |> deg2rad, site.h) #m
    # r_site_ECI = r_ecef_to_eci(PEF(), J2000(), jd) * r_site_ECEF #km
    # ρ = norm(r*1e-3 - r_site_ECI)
    ρ = norm(1e-3*r_ECEF - r_site_ECEF)
    return az, el, ρ
end

function batch_process(Az::AbstractVector, El::AbstractVector, ρ_n::AbstractVector, JD::AbstractVector, site::Site, bias:: AbstractVector; RMStol=1e-3, alg = :gibbs) # or :hgibbs
    @assert length(bias) == 3
    @assert length(Az) == length(El) && length(ρ_n) == length(Az) && length(JD) == length(Az)

    calcRMS(y, W, N, n=3) = sqrt((y' * diagm(repeat(diag(W/norm(W)), N)) * y) / (n*N))

    W = diagm(bias.^-2)
    # W = W/norm(W)
    r_site_ECEF = 1e-3*geodetic_to_ecef(site.lat |> deg2rad, site.lon |> deg2rad, site.h) #m

    ρ_NED = ρ_n .*  [ [cosd(az)*cosd(el), cosd(el)*sind(az), -sind(el)] for (az, el) in zip(Az, El) ]
    ρ_ECEF = 1e-3*ned_to_ecef.(ρ_NED*1e3, site.lat |> deg2rad, site.lon |> deg2rad, site.h)
    r_ECEF = Ref(r_site_ECEF) .+ ρ_ECEF
  
    r_ECI = r_ecef_to_eci.(PEF(), J2000(), JD) .* r_ECEF

    r_mat = genObsMat!(r_ECI)
    JD_mat = genObsMat!(JD)
    _ = genObsMat!.([Az, El, ρ_n]) # remove extra elements

    if alg == :gibbs
        v2 = [gibbs(vecs...) for vecs in eachrow(r_mat)]
    elseif alg == :hgibbs
        v2 = [herrick_gibbs(vecs..., jds) for (vecs, jds) in zip(eachrow(r_mat), eachrow(JD_mat))]
    else
        @error "Invalid algorithm, try :gibbs or :hgibbs"
    end
    SVs = OrbitStateVector.(JD_mat[:,2], 1e3*r_mat[:,2], 1e3*v2)
    r, v = unzip(Propagators.propagate_to_epoch.(Val(:TwoBody), JD[1], sv_to_kepler.(SVs)))

    # SV_nom_1 = OrbitStateVector(JD[1], mean(r), mean(v))
    # r, v = unzip(Propagators.propagate_to_epoch.(Val(:TwoBody), JD, Ref(sv_to_kepler(SV_nom_1))))
    # SV_nom_i = OrbitStateVector.(JD, r, v)
    # X_nom_i = 1e-3*Matrix(hcat([r'...;], [v'...;]))

    X_nom_1 = 1e-3*vcat(mean(r), mean(v))
    X_nom_i = propagate_state(X_nom_1, JD)

    obs_nom = state_to_obs.(eachrow(X_nom_i), JD, Ref(site))
    y = reduce(vcat, [[az, el, rho] - [obs...] for (az, el, rho, obs) in zip(Az, El, ρ_n, obs_nom)])

    H_long = Matrix{Float64}(undef, 3*length(Az), length(X_nom_1)) 

    RMS = calcRMS(y, W, size(X_nom_i, 1))
    P = nothing; HtWH = nothing; HtWy = nothing
    while true
        HtWH = zeros(Float64, (6,6))
        HtWy = zeros(Float64, 6)
        
        for (i, X_nom) in enumerate(eachrow(X_nom_i))
            δ = X_nom*.0001
            H = Matrix{Float64}(undef, 3, length(X_nom)) 
            for j in 1:6
                X_mod = Vector(copy(X_nom))
                if i != 1
                    SV = OrbitStateVector(JD[i], 1e3*X_mod[1:3], 1e3*X_mod[4:6])
                    r, v = Propagators.propagate_to_epoch(Val(:TwoBody), JD[1], sv_to_kepler(SV))
                    X_mod = Vector(1e-3*vcat(r, v)) # update to time 1
                    X_mod[j] += δ[j] # perturb at t1
                    SV = OrbitStateVector(JD[1], 1e3*X_mod[1:3], 1e3*X_mod[4:6])
                    r, v = Propagators.propagate_to_epoch(Val(:TwoBody), JD[i], sv_to_kepler(SV)) # back to ti
                    X_mod = Vector(1e-3*vcat(r, v))
                else
                    X_mod[j] += δ[j]
                end
                obs_mod = vcat(state_to_obs(X_mod, JD[i], site)...)
                H[1:3, j] .= (obs_mod - [obs_nom[i]...]) / δ[j]
            end
            HtWH += H' * W * H
            HtWy += H' * W * y[(i-1)*3 + 1:(i-1)*3 + 3]
            H_long[(i-1)*3 + 1:(i-1)*3 + 3, 1:6] .= H 
        end
        
        P = !isapprox(det(HtWH), 0) ? inv(HtWH) : pinv(HtWH) 
        x̂ = P * HtWy
        X_nom_1 = X_nom_i[1,:] + x̂

        X_nom_i = propagate_state(X_nom_1, JD)
        obs_nom = state_to_obs.(eachrow(X_nom_i), JD, Ref(site))
        y = reduce(vcat, [[az, el, rho] - [obs...] for (az, el, rho, obs) in zip(Az, El, ρ_n, obs_nom)])

        newRMS = calcRMS(y, W, size(X_nom_i, 1))
        println("RMS diff is $(abs(RMS - newRMS))")
        abs(RMS - newRMS) > RMStol || break
        RMS = newRMS
        println("New RMS is $RMS")
    end
    out = BatchData(X_nom_i[1,:], P, HtWH, HtWy, diagm(repeat(diag(W), length(Az))), H_long)
    return out
end

batch1 = batch_process(Az, El, ρ_n, JD, site, bias; alg=:hgibbs, RMStol=1e-3)
v = batch1.X_nom[4:6]; r = batch1.X_nom[1:3]

r_truth = [5748.78906896, 2679.87040097, 3442.77396002] 
v_truth = [4.32551211, -1.92323054, -5.72101322]

println("Err r is $(norm(r_truth - r)) km")
println("Err v is $(norm(v_truth - v)) km/s")

σ = sqrt.(diag(batch1.P))*1e3
eig = eigen(batch1.P[1:3, 1:3])
r_err = eig.values*1e3 #m 
r_dir = eig.vectors


# using GLMakie
# Earth = Sphere(Point3f(0), rₑ)
# F, scene = mesh(Earth; alpha=.25)
# # Makie.scatter!(scene, Point3f(r_site_ECEF); markersize=15)
# # Makie.scatter!(scene, Point3f.(ρ_ECEF); markersize=15) 
# Makie.scatter!(scene, Point3f.(r_site_ECI); markersize=15)  
# Makie.scatter!(scene, Point3f.(r_ECI); markersize=15) 
# # Makie.scatter!(scene, Point3f.(r_ECI2); markersize=15) 
# # F2, scene = Makie.scatter(Point3f.(ρ_NED))

## Problem 3

data_n = DataFrame(CSV.File("HW 2/prob3.txt"; header=false))

times_n = DateTime.(data_n[:,1], dateformat"yyyy mm dd  HH MM SS.s")
JD_n = date_to_jd.(times_n)

Az, El, ρ_n = data_n[:,2], data_n[:, 3], data_n[:, end]

batch2 = batch_process(Az, El, ρ_n, JD_n, site, bias; alg=:hgibbs, RMStol=1e-3)
v = batch2.X_nom[4:6]; r = batch2.X_nom[1:3]


σ_n = sqrt.(diag(batch2.P))*1e3
eig = eigen(batch2.P[1:3, 1:3])
r_err_n = eig.values*1e3 #m 
r_dir_n = eig.vectors

##
# W̃ = cat(batch1.W, batch2.W; dims=(1,2)) 

P0kn = batch1.P - batch1.P * batch2.H' * inv(inv(batch2.W) +  batch2.H*batch1.P*batch2.H') * batch2.H * batch1.P
x̂ = P0kn * (batch2.HtWy + batch1.HtWy)
X_nom_seq = batch1.X_nom + x̂
σ_seq = sqrt.(diag(P0kn))*1e3
eig = eigen(P0kn[1:3, 1:3])
r_err_seq = eig.values*1e3 #m 
r_dir_seq = eig.vectors

println("Err r is $(norm(r_truth - X_nom_seq[1:3])) km")
println("Err v is $(norm(v_truth - X_nom_seq[4:6])) km/s")


## Tables
# tblmat = Matrix{Any}(undef, 2, 3)
# tblmat[1, :] .= (batch1.X_nom[1:3], X_nom_seq[1:3])
# tblmat[2, :] .= (batch1.X_nom[3:6], X_nom_seq[3:6])
# pretty_table(tblmat)