using Statistics
using SatelliteToolbox
using SatelliteToolboxPropagators
using ReferenceFrameRotations
using Revise
include("JoshOrbitsLib2.jl")
plotlyjs()

doP = [false,true,true]

# prob1
if doP[1]

xo = range(1,10)|>Vector
yo = [1,2,2,3,4,4,5,7,9,9]

# test case
# xo = range(1,8)|>Vector
# yo = [1,1,2,3,3,4,4,6]

A = [ones(length(xo)) xo] # build A / H matrix for a linear system.
X̂ = inv(A'*A)*A'*yo

ϵ = [yoi - X̂[1]-X̂[2]*xoi for (yoi,xoi) in zip(yo,xo)] .|> abs
RMS = sqrt((1/length(xo))*sum(ϵ.^2))
σ = std(ϵ)
println("Before removing any points,\nthe RMS is $(RMS) and the coeffcients are $(X̂).\nThe greatest error is $(maximum(ϵ)).")

cutoff = 2σ
ind = ϵ.<= cutoff
println(ind)
xo = xo[ind]
yo = yo[ind]

A = [ones(length(xo)) xo] # build A / H matrix for a linear system.
X̂ = inv(A'*A)*A'*yo

ϵ = [yoi - X̂[1]-X̂[2]*xoi for (yoi,xoi) in zip(yo,xo)]
RMS = sqrt((1/length(xo))*sum(ϵ.^2))
σ = std(ϵ)
println("After removing any points that have an error greater than the previous $(cutoff),\nthe RMS is $(RMS) and the coeffcients are $(X̂).\nThe greatest error is $(maximum(ϵ)).")

end
# Q2
if doP[2]

d = [
1995 01 29 02 38 37.00 60.4991 16.1932 2047.50200
1995 01 29 02 38 49.00 62.1435 17.2761 1984.67700
1995 01 29 02 39 02.00 64.0566 18.5515 1918.48900
1995 01 29 02 39 14.00 65.8882 19.7261 1859.32000
1995 01 29 02 39 26.00 67.9320 20.9351 1802.18600
1995 01 29 02 39 38.00 70.1187 22.1319 1747.29000
1995 01 29 02 39 50.00 72.5159 23.3891 1694.89100
1995 01 29 02 40 03.00 75.3066 24.7484 1641.20100
1995 01 29 02 40 15.00 78.1000 25.9799 1594.77000
1995 01 29 02 40 27.00 81.1197 27.1896 1551.64000
]

#                   biases↓                 noise↓
#                   ρ    
#Kaena Point, Ha    80.0 0.0081 0.0045      92.5 0.0224 0.0139

jd = [DateTime(i[1:6]...) for i in eachrow(d)]
jd = date_to_jd.(jd)
az = d[:,7] .|> deg2rad
el = d[:,8] .|> deg2rad
ρ = d[:,9]


Φ = 21.5721  |> deg2rad   # latitude
λ = -158.2666  |> deg2rad # longitude
h = 300.2/1000            # altitude (m)

# ρb = rhoElAzToECI.(ρ.-80/1000,el.-deg2rad(0.0045),az.-deg2rad(0.0081),Φ,λ,h,jd)
ρv =  rhoElAzToECI.(ρ,el,az,Φ,λ,h,jd) 


Rs = [Φ,λ,h]
RsECEF = geodetic_to_ecef(Rs...) # take R site to ECEF
q=r_ecef_to_eci.(Quaternion,PEF(),J2000(),jd)
RsECI = vecRotate.(Ref(RsECEF),q)./1000 # rotate into ECI frame
r = RsECI+ρv .|> Vector # list of ten rs as found in each obs
# rb = RsECI+ρb

rss = [[r[1],r[3],r[7]],[r[2],r[4],r[8]],[r[3],r[5],r[9]]] # 9 rs (skipped #6) grouped into 3s for IOD
jdss = [[jd[1],jd[3],jd[7]],[jd[2],jd[4],jd[8]],[jd[3],jd[5],jd[9]]] # the associated times ↑
v_hg = doHgibbs.(rss,jdss) 
v_hg = scrape.(v_hg) # the vs associated with the middle of each set ↑
r_hg = [r[3],r[4],r[5]] # middle of each set | input to statistics part
# println(jd)
jd_hg = [jd[3],jd[4],jd[5]]
# println(jd)

# el = [el[3],el[4],el[5]]
# az = [az[3],az[4],az[5]]
# ρ = [ρ[3],ρ[4],ρ[5]]
# RsECI = [RsECI[3],RsECI[4],RsECI[5]]


orb₀ = rv_to_kepler.(r_hg*1000,v_hg*1000,jd_hg)
# println(kepler_to_rv(orb₀[1])./1000)
props = Propagators.init.(Val(:TwoBody),orb₀)
rv = Propagators.propagate_to_epoch!.(props,jd[1])
r_hg = scrape.(rv,1)./1000
v_hg = scrape.(rv,2)./1000

function propagateX̂(X̂,jd0,jd)
    r = X̂[1:3]
    v = X̂[4:6]
    orb₀ = rv_to_kepler(r*1000,v*1000,jd0)
    props = Propagators.init(Val(:TwoBody),orb₀)
    rv = Propagators.propagate_to_epoch!(props,jd)
    r = rv[1]./1000
    v = rv[2]./1000
    X̂ = vcat(r,v) |> Vector
end

function predictObs(X̂,RsECI,jd) # grab the constants from the workspace, its fine
    r = X̂[1:3]
    elc,azc = calcElAz(Φ,λ,r-RsECI,jd)
    ρc = norm(r-RsECI)
    return elc,azc,ρc
end


# nom vectors going into the iteration
r_hg = [sum(scrape.(r_hg,1)),sum(scrape.(r_hg,2)),sum(scrape.(r_hg,3))]./length(r_hg)
v_hg = [sum(scrape.(v_hg,1)),sum(scrape.(v_hg,2)),sum(scrape.(v_hg,3))]./length(v_hg)

X̂ₙ = vcat(r_hg,v_hg)
# println(X̂ₙ)


δρ,δaz,δel = [92.5/1000,deg2rad(0.0224),deg2rad(0.0139)]
W = 
[1/(δρ^2) 0 0
0 1/(δaz^2) 0
0 0 1/(δel^2)]
# println(W)

numMeas = 3 # dimension of observations
N = length(el) # numeber of observations
n = 6 # dimension of state vector

δratio = .0001 # amount we mod by 
rmsNew = Inf
iterations = 0
X̂out = []
rmslong = []
while true
    global iterations
    global rmsNew
    global rmsOld
    global X̂ₙᵢ 
    global X̂ₙ 
    global b̃ = []
    global X̂out
    global AWA 
    global AWb 
    global rmslong
    AWA = zeros(n,n)
    AWb = zeros(n)

    for i in 1:N # loop over the 10 observations   
        X̂ₙᵢ = propagateX̂(X̂ₙ,jd[1],jd[i])
        elc,azc,ρc = predictObs(X̂ₙᵢ,RsECI[i],jd[i]) # based on our best guess, what do we think the observations should be
        r̄ = [ρ[i]-ρc,az[i]-azc,el[i]-elc] # calc residuals based on this, its a measure of how bad our guess is
        # println(r̄)
        append!(b̃,r̄)

        # build the A matrix for ti
        A = zeros(3,6)
        for j in 1:n # loop over the 6 elements of SV
            X̂ₘ = copy(X̂ₙ)
            δ = X̂ₙ[j] * δratio
            X̂ₘ[j] = X̂ₙ[j] + δ # take our best guess and muck it up a bit

            X̂ₘ = propagateX̂(X̂ₘ,jd[1],jd[i]) # see where we would be at ti if we did this        
            elₘ,azₘ,ρₘ = predictObs(X̂ₘ,RsECI[i],jd[i]) # for that vecotr above, what would the observations be in a perfect world
            elₙ,azₙ,ρₙ = predictObs(X̂ₙᵢ,RsECI[i],jd[i]) # see what the nominal would look like at the same time
            A[1,j] = ( elₘ - elₙ ) / δ # build A 
            A[2,j] = ( elₘ - azₙ ) / δ
            A[3,j] = ( ρₘ - ρₙ ) / δ
        end
        AWA += A'*W*A # HtWH
        AWb += A'*W*r̄ # HtWy
    end
    try
        global P = inv(AWA)
    catch
        global P = pinv(AWA)
    end

    Wlong = diagm(repeat(diag(W),N))

    x̂ = P*AWb
    X̂ₙ += x̂
    rmsOld = copy(rmsNew)
    rmsNew = sqrt((b̃'*Wlong*b̃)/(numMeas*N*norm(W)))

    iterations+=1    
    
    println()
    println(iterations)
    # println(x̂)
    println(rmsNew)
    println(X̂ₙ)
    # println(AWb)
    # println(AWA)
    push!(X̂out,X̂ₙ)
    append!(rmslong,rmsNew)
    
    # abs(rmsOld-rmsNew)/rmsNew > .0001 || break
    iterations < 200 || break
end


plot()
# plot!(1:length(X̂out)|>Vector,scrape.(X̂out,1))
# plot!(1:length(X̂out)|>Vector,scrape.(X̂out,2))
# plot!(1:length(X̂out)|>Vector,scrape.(X̂out,3),show = true, legend = false)
# plot(1:length(rmslong)|>Vector,rmslong,show = true, legend = false)




# while true #RMSprev-RMScur/RMSprev < tol (1e-3)
#     # b̃         - residual matrix
#     # A         - partial derivative matrix
#     # W         - weights matrix
#     # δx̂        - change in x̂
#     # P = A'WA  - covariance matrix

#     # Propagate the nominal state to each ti and find computed observations
#     for i in enumerate(r)
        
#         elc,azc = calcElAz(Φ,λ,r-RsECI,jd)
#         ρc = norm(r-RsECI)
#         r̄ = [ρ-ρc,az-azc,el-elc]

#         # calc weights
#         W = 
#         [[1/δρ^2,0,0],
#         [0,1/δaz^2,0],
#         [0,0,1/δel^2]]
#     end
# end

end
# Q3
if doP[3]
# 1995 01 29 02 40 39.00 84.3708 28.3560 1512.08500
# 1995 01 29 02 40 51.00 87.8618 29.4884 1476.41500
# 1995 01 29 02 41 03.00 91.5955 30.5167 1444.91500
# 1995 01 29 02 41 15.00 95.5524 31.4474 1417.88000
# 1995 01 29 02 41 27.00 99.7329 32.2425 1395.56300
# 1995 01 29 02 41 39.00 104.0882 32.8791 1378.20200
# 1995 01 29 02 41 51.00 108.6635 33.3788 1366.01000
# 1995 01 29 02 42 03.00 113.2254 33.5998 1359.10000





end
println("End of line")