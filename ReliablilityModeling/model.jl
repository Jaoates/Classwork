using Plots
using Revise
plotlyjs()

function R(α,β,t) exp(-α*t)^β end

###############################
λ = 1/6
function Relectrical(t) R(λ,1,t) end

λ = 1/6
function Rlube(t) R(1/6,1,t) end

α = 1/6
β = 1.2
function Rmotor(t) R(1/5,1.2,t) end

α = 1/6
β = 1.2
function Rbearings(t) R(1/5,1.2,t) end

α = 1/6
β = 1.2
function Rspline(t) R(1/5,1.2,t) end
###############################

function Rₛ(t) prod([Relectrical(t),Rlube(t),Rmotor(t),Rbearings(t),Rspline(t)]) end

t = range(0,10,100)
plot()
for i in range(0,2,20)
    plot!(t,R.(.5,i,t))
end
plot!(show = true)