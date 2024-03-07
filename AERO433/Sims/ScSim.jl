using Revise
includet("SimLib.jl")
using .PowerFunction
using PlotlyJS

r = [6539, 2708,0]
v = [0.4,-.96,7.4]
X = vcat(r,v)
jd = 0:.01:1
params = SimParams()
# p = Psar(rk(r,v,jd,params)...)
# p = Psun(rk(r,v,jd,params)...)
# p = Pcom(rk(r,v,jd,params)...)
# p = Pheater(rk(r,v,jd,params)...)
# p = Pobc(rk(r,v,jd,params)...)

p = calcPower.(jd,Ref(X),Ref(params))
