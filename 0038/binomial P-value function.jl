# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10)

function pvalue_bin_wilson(k, p; n=20)
    z = (k - n*p)/√(n*p*(1 - p))
    2ccdf(Normal(), abs(z))
end

# %%
n = 20
α = 0.05
x = range(-0.5, n+0.5, 500)
p = range(0, 1, 500)
heatmap(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n) ≥ α;
    clim=(0, 1), colorbar=false)
plot!(xguide="data k", yguide="parameter p")
plot!(xtick=0:n, ytick=0:0.1:1)
title!("Wilson's confidence intervals for n = $n, α = $α")
plot!(size=(580, 500))

# %%
n = 20
x = range(-0.5, n+0.5, 500)
p = range(0, 1, 500)
heatmap(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); clim=(0, 1))
plot!(xguide="data k", yguide="parameter p")
plot!(xtick=0:n, ytick=0:0.1:1)
title!("Wilson's P-value function for n = $n")
plot!(size=(640, 500))

# %%
n, k = 20, 6
plot(p, p -> pvalue_bin_wilson(k, p; n); label="")
plot!(xguide="parameter p", yguide="P-value")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("Wilson's P-value function for n=$n, k=$k")

# %%
n = 20
x = range(-0.5, n+0.5, 211)
p = range(0, 1, 201)
@time anim = @animate for t in 0:3:359
    surface(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); colorbar=false, camera=(t, 60))
    plot!(xguide="k", yguide="p", zguide="P-value")
    plot!(size=(640, 600), zlim=(0, 1))
end
gif(anim, "Wilson's P-value function n=$n (60).gif")

# %%
n = 20
x = range(-0.5, n+0.5, 211)
p = range(0, 1, 201)
@time anim = @animate for t in 0:3:359
    surface(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); colorbar=false, camera=(t, 45))
    plot!(xguide="k", yguide="p", zguide="P-value")
    plot!(size=(640, 600), zlim=(0, 1))
end
gif(anim, "Wilson's P-value function n=$n (45).gif")

# %%
