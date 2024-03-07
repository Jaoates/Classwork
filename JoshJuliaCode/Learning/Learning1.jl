using LinearAlgebra
using Plots
using DifferentialEquations
using StaticArrays

function hello(myName)
    println("hi there $myName")
    print()
end

# initialize vector
v1 = [1,2,3]

# transpose vector
v2 = v1'

# vector multiplication
M1 = v1*v2

# element wise multiplication
v3 = v1.*v1

# cross product
cross(v1,v3)

f(x)=sin.(x)
g(x)=cos.(x)

x = 0:.1:10
println(typeof(x))
# x = collect(x)
y1 = f(x)
y2 = g(x)

plot(x,y1)
plot!(x,y2)

scatter(y2,y1,x,xaxis = "x",yaxis = "y",zaxis = "z")

bar(y2)

vout = [0,0,0]

for (i,v) in enumerate(v1)
    vout[1] = 2*v
end

gr()

#Half-life of Carbon-14 is 5,730 years.
C₁ = 5.730

#Parameters
ω = .5

# #Initial Conditions
# x₀ = [1.0]
# dx₀ = [π / 2]
# tspan = (0.0, 100)

# ϕ = atan((dx₀[1] / ω) / x₀[1])
# A = √(x₀[1]^2 + dx₀[1]^2)


# #Define the problem
# function myFun(ddu, du, u, ω, t)
#     ddu .= -ω^2 * u
# end
# # Pass to solver
# prob = SecondOrderODEProblem(myFun, dx₀, x₀, tspan, ω)
# sol = solve(prob)
# plot(sol)


#Initial Conditions
x₀ = 1.0
dx₀ = 1.5
tspan = (0.0, 2)

#Define the problem
X₀ = [x₀,dx₀] 

# X₀ = SVector(X₀,2)
function myFun!(dX, X, ω, t)
    print(X)
    v = X[1]
    dv = X[2]
    ddv = -ω^2 * v
    dX = [dv,ddv]
    # print(t)
    println(dX)
    return dX
end

#Pass to solver
prob = ODEProblem(myFun!, X₀, tspan, ω)
sol = solve(prob, saveat=.5, reltol=1e-8)

#Plot
plot(sol, linewidth = 2)
# plot!(sol.t, t -> exp(-C₁ * t), lw = 3, ls = :dash, label = "Analytical Solution")

# function lorenz!(du, u, p, t)
#     du[1] = 10.0 * (u[2] - u[1])
#     du[2] = u[1] * (28.0 - u[3]) - u[2]
#     du[3] = u[1] * u[2] - (8 / 3) * u[3]
# end

# using DifferentialEquations
# u0 = [1.0; 0.0; 0.0]
# tspan = (0.0, 100.0)
# prob = ODEProblem(lorenz!, u0, tspan)
# sol = solve(prob)