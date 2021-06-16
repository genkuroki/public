using Pkg; Pkg.activate("HeatEqCalc")
using HeatEqCalc

ENV["GKSwstype"]=100
using Plots
using Printf

sol = HeatEqCalc.load_sol("sol")
t, x, u, m = sol
u0 = u[:, begin]

P = plot(; legend=false)
for k in 1:m
    plot!(x, (p -> p[k]).(u0); lw=0.1, color=:blue)
end
plot()
savefig(P, "heateq0.png")

ylim = (minimum(minimum.(u0)), maximum(maximum.(u0)))
anim = @animate for i in 1:length(t)
    title = @sprintf("t = %4.2f", t[i])
    plot(; legend=false)
    for k in 1:m
        plot!(x, (p -> p[k]).(u[:, i]); ylim, lw=0.1, color=:blue)
    end
    title!(title)
end
gif(anim, "heateq.gif")
