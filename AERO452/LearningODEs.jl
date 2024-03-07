using Plots
using ReferenceFrameRotations
using DifferentialEquations
using Revise

plotlyjs()

function EOM(X,p,t)
    θ = X[1]
    ω = X[2]
    g,k = p
    α = -sin(θ)*g -ω*k
    Ẋ= [ω,α]
end

p = [9.8,1]
X₀ = [π/2,0]
tspan = (0,10)
myProb = ODEProblem(EOM,X₀,tspan,p)

sol = solve(myProb,)

tspan = 0:.01:5
solt = sol(tspan)

u = [solt[1,:] solt[2,:]]

plot(eachcol(u)...,title = "State Space",show = true)




function EOM!(Ẋ,X,p,t)
    θ = X[1]
    ω = X[2]
    g,k = p
    α = -sin(θ)*g -ω*k
    Ẋ[1] = ω
    Ẋ[2] = α
    nothing
end

p = [9.8,1]
X₀ = [π/2,0]
tspan = (0,10)
myProb = ODEProblem(EOM!,X₀,tspan,p)
sol = solve(myProb)

tspan = 0:.01:5
solt = sol(tspan)

u = [solt[1,:] solt[2,:]]

plot(eachcol(u)...,title = "State Space Again")