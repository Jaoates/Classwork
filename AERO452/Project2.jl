# second project - Josh Oates and Reno Brown

println("Program Start")

include("JoshOrbitsLib.jl")
include("OtherOrbitsLib.jl")
using SatelliteToolbox
using Plots
using ReferenceFrameRotations
using DifferentialEquations
using Revise
using Unzip
using FileIO
using JLD2

## Settings
plotlyjs()
doPlots = true
##

## Initialize Things
satelliteNames = ["GOES 18","GRACE-FO 1","NOAA 21"]
# fetching code 
    # fetcher = create_tle_fetcher(CelestrakTleFetcher)
    # tles = map(x-> fetch_tles(fetcher,satellite_name = x)[1],satelliteNames)
    # FileIO.save("xtlesx.jld2","tles",tles)
tles = FileIO.load("tles.jld2")["tles"]
props = map(x->Propagators.init(Val(:SGP4),x),tles)
jd₀= date_to_jd(2023,12,3)
rv₀ = Propagators.propagate_to_epoch!.(props,jd₀)
r₀ = map(x->x[1]|>Vector,rv₀)
v₀ = map(x->x[2]|>Vector,rv₀)
orb₀ = rv_to_kepler.(r₀,v₀)
propsOSC = map(x->Propagators.init(Val(:TwoBody), x),orb₀)



r₀/=1000
v₀/=1000
X₀ = map((x,y)->[x;y],r₀,v₀ )

struct Params
    μn
    jd₀
    BC
    SC
end
Params(BC,SC) = Params(μsun,jd₀,BC,SC)

m = [5192,487,2930]
A = [34.16,6.064866,22.7612448]

Cr = 2
SC = m./(Cr.*A)

Cd = 2.2
BC = m./(Cd.*A)

params = Params.(BC,SC)

function EOM1(X,p,t)
    μsun = p.μn
    jd₀ = p.jd₀
    BC = p.BC
    SC = p.SC
    r = X[1:3]
    v = X[4:6]
    jd = jd₀ + t/(24*3600)

    rsun = sun_position_mod(jd)/1000
    Q = r_eci_to_eci.(Quaternion,MOD(),J2000(),jd)
    rsun = vecRotate(rsun,Q)

    a = calcAccel(r) + calcNBodyAccel(r,rsun,μsun) + sum(calcJaccel(r)) + calcSRPaccel(r,rsun,SC) + calcDrag(r,v,BC) 
    Ẋ = [v;a]
end

function EOM2(X,p,t)
    μsun = p.μn
    jd₀ = p.jd₀
    #BC = p.BC
    SC = p.SC
    r = X[1:3]
    v = X[4:6]
    jd = jd₀ + t/(24*3600)

    rsun = sun_position_mod(jd)/1000
    Q = r_eci_to_eci.(Quaternion,MOD(),J2000(),jd)
    rsun = vecRotate(rsun,Q)

    a = calcAccel(r) + calcNBodyAccel(r,rsun,μsun) + sum(calcJaccel(r)) + calcSRPaccel(r,rsun,SC)# + calcDrag(r,v,BC) 
    Ẋ = [v;a]
end

EOM = [EOM2,EOM1,EOM1]
##

## Uncorrected Propagation
tspan = [1,1,1]*3600*24*(365/4)

prob = ODEProblem.(EOM,X₀,tspan,params)
sol =  @time solve.(prob,reltol = 1e-10,abstol = 1e-10)
r = map(x->x[1:3, :],sol)
v = map(x->x[4:6, :],sol)
t = map(x->x.t,sol)
r = map(x->collect(eachcol(x)),r)
v = map(x->collect(eachcol(x)),v)

orbs = map((r,v)->rv_to_kepler.(r.*1000,v.*1000,0),r,v)

function plotRelvantInfo(orb,r,t,title)

    # t/=(24*3600)
    plot(scrape.(r,1),scrape.(r,2),scrape.(r,3),show = false)
    plotEarth!(false)
    plot!(title = "trajectory over time span $(title)",xlabel = "x",ylabel = "y",zlabel = "z",show = doPlots)
    
    u = vcat(Vector.(CoeStructure.(orb)))
    plot(t,rad2deg.(scrape.(u,3).-u[1][3]),title="ΔInc vs t $(title)",ylabel="Deg",xlabel="s",show=doPlots)
    plot(t,rad2deg.(scrape.(u,4).-u[1][4]),title="ΔΩ vs t $(title)",ylabel="Deg",xlabel="s",show=doPlots)
    plot(t,rad2deg.(scrape.(u,5).-u[1][5]),title="Δω vs t $(title)",ylabel="Deg",xlabel="s",show=doPlots)
    # plot(t,rad2deg.(scrape.(u,6).-u[1][6]),title="Anomally vs t $(title)",ylabel="Deg",xlabel="s",show=true)
    RaRp = calcRaRp.(orb)
    plot(t,scrape.(RaRp,1).-re,label="za",show=doPlots)
    plot!(t,scrape.(RaRp,2).-re,label="zp",title="za and zp vs t $(title)",ylabel="km",xlabel="s",show=doPlots)
end

# plotRelvantInfo.(orbs,r,t,satelliteNames)
##

## Corrections Setup
function plotRelvantInfo2(orb₀,r,v,t,title)
    orb = rv_to_kepler.(r.*1000,v.*1000,0)
    # t/=(24*3600)
    plot(scrape.(r,1),scrape.(r,2),scrape.(r,3),show = false)
    plotEarth!(false)
    plot!(title = "trajectory over time span $(title)",xlabel = "x",ylabel = "y",zlabel = "z",show = doPlots)
    
    u = vcat(Vector.(CoeStructure.(orb)))
    plot(t,rad2deg.(scrape.(u,3).-u[1][3]),title="ΔInc vs t $(title)",ylabel="Deg",xlabel="s",show=false)
    # plot!([t[1],t[end]],[1,1]*rad2deg(orb₀.i),show=doPlots)
    plot(t,rad2deg.(scrape.(u,4).-u[1][4]),title="ΔΩ vs t $(title)",ylabel="Deg",xlabel="s",show=doPlots)
    # plot!([t[1],t[end]],[1,1]*rad2deg(orb₀.Ω),show=doPlots)
    plot(t,rad2deg.(scrape.(u,5).-u[1][5]),title="Δω vs t $(title)",ylabel="Deg",xlabel="s",show=doPlots)
    # plot!([t[1],t[end]],[1,1]*rad2deg(orb₀.ω),show=doPlots)

    # plot(t,rad2deg.(scrape.(u,6).-u[1][6]),title="Anomally vs t $(title)",ylabel="Deg",xlabel="s",show=true)
    
    RaRp = calcRaRp(orb₀)
    plot([t[1],t[end]],[1,1]*(RaRp[1].-re),label="za osc",show=false)
    plot([t[1],t[end]],[1,1]*(RaRp[2].-re),label="za osc",show=false)
    RaRp = calcRaRp.(orb)
    plot(t,scrape.(RaRp,1).-re,label="za",show=false)
    plot!(t,scrape.(RaRp,2).-re,label="zp",title="za and zp vs t $(title)",ylabel="km",xlabel="s",show=doPlots)
end

function doSim(EOM,X₀,p,tEndSim,Tl,Tc,orb₀,satelliteName)
    # println("--$(satelliteName)--")
    X = [X₀]
    # p = params[i]
    # EOMi = EOM[i]
    # orbi₀ = orb₀[i]
    t = [0.0]
    # Tl = 300.0
    # Tc = (24/4)*3600
    ΔV = 0.0
    # tEndSim = 24*3600
    tol = 1e-11
    saveTime = 1
    # Corrections sim
    # println("dear god")
    i = 0
    while t[end]<tEndSim
        i+=1
        # println("while?$(i)")
        local v
        local r
        local prob
        local sol
        
        # global X

        # propoagate the perturbed orbit
        prob = ODEProblem(EOM,X[end],(t[end],t[end]+Tc),p)
        sol =  solve(prob,reltol = tol,abstol = tol,saveat = saveTime)
        
        # println(t[end])
        # println(X[end])
        # println(sol.u[1])
        # len = length(X)

        append!(X,sol.u)
        # println.(X[len:len+5])
        append!(t,sol.t)
        # println(t[end])

        # create an osculating orbit and find out where you want to be at the end of the transfer
        r,v = Propagators.propagate(Val(:TwoBody),t[end]+Tl,orb₀)
        r/=1000
        v/=1000

        # solve lamberts to get your dvs and how to leave current orbit and arrive at osculating
        v1,v2=lambert_reno(X[end][1:3],r,Tl)
        ΔV+=norm(X[end][4:6] - v1)
        # X[end][4:6] = v1
        
        # println(t[end])
        push!(X,vcat(X[end][1:3],v1))
        push!(t,t[end])
        # println(t[end])

        # propogate the perturbed orbit after starting transfer
        prob = ODEProblem(EOM,X[end],(t[end],t[end]+Tl),p)
        sol =  solve(prob,reltol = tol,abstol = tol,saveat = saveTime)
        
        # println(t[end])
        # println(X[end])
        # println(sol.u[1])
        # len = length(X)

        append!(X,sol.u)
        # println.(X[len:len+5])
        append!(t,sol.t)

        ΔV+=norm(X[end][4:6] - v)

        # burn at end of transfer to resume osculating orbit
        # X[end][4:6] = v
        # println(t[end])
        push!(X,vcat(X[end][1:3],v))
        push!(t,t[end])
        # println(t[end])

    end

    # if length(t)<length(X)
    #     X = X[1:length(t)]
    # elseif length(t)>length(X)
    #     t = t[1:length(X)]
    # end

    # r = map(x->x[1:3],X)
    # v = map(x->x[4:6],X)
    # orbs = rv_to_kepler.(r.*1000,v.*1000,0)
    # plotRelvantInfo2(orb₀,r,v,t,satelliteName)

    # println("ΔV: $(ΔV*1000)")
    ΔV    
end

## finding Min ΔVs

#Find the min ΔV for GRACE
    Tl = range(200,1200,50)
    Tc = range(10,600,50)

    i = range(1,length(Tl))
    j = range(1,length(Tc))
    ΔV = zeros(length(Tl),length(Tc))

    for i in i
        for j in j
            println("i = $(i)")
            println("j = $(j)")
            tspan = [1,1,1]*3600*24*(2)
            ΔV[i,j] = doSim(EOM[2],X₀[2],params[2],tspan[2],Tl[i],Tc[j],orb₀[2],satelliteNames[2])
            println("ΔV: $(ΔV[i,j]*1000)")
        end
    end
    m,i=findmin(ΔV)
    FileIO.save("mindvGRACE.jld2","mindv",(m,i,ΔV,Tl,Tc))


#Find the min ΔV for GOES-18
    Tl = range(200,1200,50)
    Tc = range(10,600,50)

    i = range(1,length(Tl))
    j = range(1,length(Tc))
    ΔV = zeros(length(Tl),length(Tc))

    for i in i
        for j in j
            println("i = $(i)")
            println("j = $(j)")
            tspan = [1,1,1]*3600*24*(2)
            ΔV[i,j] = doSim(EOM[1],X₀[1],params[1],tspan[1],Tl[i],Tc[j],orb₀[1],satelliteNames[1])
            println("ΔV: $(ΔV[i,j]*1000)")
        end
    end
    m,i=findmin(ΔV)
    FileIO.save("mindvGOES18.jld2","mindv",(m,i,ΔV,Tl,Tc))

m,i,ΔV,Tl,Tc = FileIO.load("mindvGOES18.jld2")["mindv"]
timePairs = [[xi, yi] for xi in Tl, yi in Tc]
Tl = Tl[i[1]]
Tc = Tc[i[2]]
surface(scrape.(timePairs,1),scrape.(timePairs,2),ΔV,show = true,xlabel = "Tl [s]",ylabel = "Tc [s]",zlabel = "ΔV [km/s]",title = "GOES-18 minimum ΔV")
scatter3d!([Tl],[Tc],[m],label = "Minimum")
contour(scrape.(timePairs,1)|>vec,scrape.(timePairs,2)|>vec,ΔV|>vec,levels=20,fill = true,show = true,xlabel = "Tl [s]",ylabel = "Tc [s]",title = "GOES-18 minimum ΔV")
scatter!([Tl],[Tc],label = "Minimum")
println("--GOES18--")
println("The minimum ΔV is $(m)km with Tl = $(Tl)s and Tc = $(Tc)s")
# hope it works
m,i,ΔV,Tl,Tc = FileIO.load("mindvGRACE.jld2")["mindv"]
timePairs = [[xi, yi] for xi in Tl, yi in Tc]
Tl = Tl[i[1]]
Tc = Tc[i[2]]
surface(scrape.(timePairs,1),scrape.(timePairs,2),ΔV,show = true,xlabel = "Tl [s]",ylabel = "Tc [s]",zlabel = "ΔV [km/s]",title = "GRACE-FO 1 minimum ΔV")
scatter3d!([Tl],[Tc],[m],label = "Minimum")
contour(scrape.(timePairs,1)|>vec,scrape.(timePairs,2)|>vec,ΔV|>vec,levels=20,fill = true,show = true,xlabel = "Tl [s]",ylabel = "Tc [s]",title = "GRACE-FO 1 minimum ΔV")
scatter!([Tl],[Tc],label = "Minimum")
println("--GRACE-FO 1--")
println("The minimum ΔV is $(m)km with Tl = $(Tl)s and Tc = $(Tc)s")


## Archive
    ## Corrections 
    # function doStep(X,params,t,tstep,EOM)
    #     println("do step")
    #     X₀ = X[end]
    #     Xold = X
    #     told = t
    #     tspan = (t[end],tstep+t[end])
    #     prob = ODEProblem(EOM,X₀,tspan,params)
    #     sol =  solve(prob,reltol = 1e-10,abstol = 1e-10)
    #     Xnew = sol.u
    #     tnew = sol.t
    #     println("told: $(told[1]),$(told[end])")
    #     println("tnew: $(tnew[1]),$(tnew[end])")
    #     # println(append!(told,tnew))
    #     append!(Xold,Xnew)
    #     append!(told,tnew)
    #     println("append: $(told[1]),$(told[end])")
    #     return Xold,told
    # end

    # function doCorrection(prop,X,params,t,tstep,ΔV,EOM)
    #     println("do correction")
    #     X₀ = X[end]

    #     r₀ = X₀[1:3]
    #     v₀ = X₀[4:6]
    #     poop = prop
    #     rf,vf = Propagators.propagate!(poop,t[end]+tstep)

    #     rf/=1000
    #     vf/=1000

    #     v1,v2=lambert_reno(r₀,rf,tstep)
    #     # println(X[end])
    #     X[end][4:6] = v1
    #     # println(X[end])
    #     X,t = doStep(X,params,t,tstep,EOM)
    #     v2 = X[end][4:6]
        
    #     ΔV += norm(v2-vf)+norm(v1-v₀)
    #     # append!(Xold,X)
    #     # append!(told,t)
    #     return X,t,ΔV
    # end

    # tstepDrift = [1,1,1]*(3000) |>Vector{Float64}
    # tstepCorrect = [1,1,1]*(3000) |>Vector{Float64}
    # tEndSim = [1,1,1]*40000 |>Vector{Float64}

    # function doSim(propsOSC,X₀,tstepDrift,tstepCorrect,tEndSim,params,EOM)
    #     ΔV = 0.0
    #     t = [0.0]
    #     X = [X₀]
    #     # println(params)
    #     while t[end]<tEndSim
    #         X,t = doStep(X,params,t,tstepDrift,EOM)
            
    #         r = map(x->x[1:3],X)
    #         v = map(x->x[4:6],X)
    #         orbs = rv_to_kepler.(r.*1000,v.*1000,0)
    #         plotRelvantInfo(orbs,r,t,"poo1")

    #         X,t,ΔV = doCorrection(propsOSC,X,params,t,tstepCorrect,ΔV,EOM)

    #         r = map(x->x[1:3],X)
    #         v = map(x->x[4:6],X)
    #         orbs = rv_to_kepler.(r.*1000,v.*1000,0)
    #         plotRelvantInfo(orbs,r,t,"poo2")
    #         # println("time: $(t[end]/24/3600)")
    #         # println("size X: $(sizeof(X))")
    #     end
    #     return X,t,ΔV
    # end

    # X,t,ΔV = doSim(propsOSC[1],X₀[1],tstepDrift[1],tstepCorrect[1],tEndSim[1],params[1],EOM[1])


    ##

    # SRP test 
        # jd₀ = 2438400.5
        # h = 63383.4
        # ecc = .025422
        # Ω = 45.3812 |> deg2rad
        # inc = 88.3924|> deg2rad
        # ω = 227.493|> deg2rad
        # θ = 343.427|> deg2rad

        # a = calcSemimajor(h,ecc,input = "he")
        # orb₀ = KeplerianElements(0,a*1000,ecc,inc,Ω,ω,θ)
        # r₀,v₀ = kepler_to_rv(orb₀)
        # r₀ /=1000 
        # v₀ /=1000 
        # X₀ = [r₀;v₀] |> Vector

        # Cr = 2
        # A_m = 2 #m^2/kg
        # SC = 1/(Cr*A_m)

        # struct Params
        #     SC
        #     jd₀
        # end
        # p = Params(SC,jd₀)

        # function EOM(X,p,t)
        #     jd₀ = p.jd₀
        #     SC = p.SC
        #     r = X[1:3]
        #     v = X[4:6]
        #     jd = jd₀ + t/(24*3600)

        #     rsun = sun_position_mod(jd)/1000
        #     Q = r_eci_to_eci.(Quaternion,MOD(),J2000(),jd)
        #     rsun = vecRotate(rsun,Q)

        #     a = calcAccel(r)+calcSRPaccel(r,rsun,SC)
        #     # println(jd-jd₀)
        #     # println(norm(calcSRPaccel(r,rsun,SC)))
        #     Ẋ = [v;a]
        # end

        # tspan = (0,1000*24*3600)
        # myProb = ODEProblem(EOM,X₀,tspan,p)
        # sol = solve(myProb,reltol = 1e-10,abstol = 1e-10)
        # r = [sol[1,:] sol[2,:] sol[3,:]]
        # v = [sol[4,:] sol[5,:] sol[6,:]]
        # t = sol.t
        # t /= (24*3600)
        # tf = t[end]
        # orbs = rv_to_kepler.(eachrow(r.*1000),eachrow(v.*1000),0)
        # u = vcat(Vector.(CoeStructure.(orbs)))

        # plot(t,rad2deg.(scrape.(u,3).-u[1][3]),title="ΔInc vs t (SRP)",ylabel="Deg",xlabel="Days",show=true)
        # plot(t,rad2deg.(scrape.(u,4).-u[1][4]),title="ΔΩ vs t (SRP)",ylabel="Deg",xlabel="Days",show=true)
        # plot(t,rad2deg.(scrape.(u,5).-u[1][5]),title="Δω vs t (SRP)",ylabel="Deg",xlabel="Days",show=true)
        # plot(eachcol(r)...,show = false)
        # plotEarth!(true)
    ##

println("Program End")