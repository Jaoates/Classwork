using Statistics
using SatelliteToolbox
using SatelliteToolboxPropagators
using ReferenceFrameRotations
using Revise
using LinearAlgebra
include("JoshOrbitsLib2.jl")


#60.4991 16.1932 2047.50200


az = 60.4991 .|> deg2rad
el = 16.1932 .|> deg2rad
ρ = 2047.50200

jd = date_to_jd(DateTime(1995,01,29,02,38,37.00))

Φ = 21.5721  |> deg2rad   # latitude
λ = -158.2666  |> deg2rad # longitude
h = 300.2/1000          # altitude (m)
h = 0

# ρb = rhoElAzToECI.(ρ.-80/1000,el.-deg2rad(0.0045),az.-deg2rad(0.0081),Φ,λ,h,jd)
ρ = rhoElAzToECI(ρ,el,az,Φ,λ,h,jd) 


Rs = [Φ,λ,h]
RsECEF = geodetic_to_ecef(Rs...) # take R site to ECEF
q=r_ecef_to_eci(Quaternion,PEF(),J2000(),jd)
RsECI = vecRotate(RsECEF,q)./1000 # rotate into ECI frame
r = RsECI+ρ |> Vector # list of ten rs as found in each obs

println(r)

el,az = calcElAz(Φ,λ,r-RsECI,jd) .|>rad2deg

lst = calcLst(λ,jd)
# el,az = rAZEL(r-RsECI,Φ,λ,lst).|>rad2deg # this function is coded up and ready to use if you need


println("$(az) | $(el)")
