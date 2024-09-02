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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, tickfontsize=6, titlefontsize=12, legendfontsize=12)

function plot_paths(paths, t, p; kwargs...)
    n = length(paths[1])
    s = n ÷ 1000
    ns = [1; s:s:n-s; n]
    plot()
    for i in 1:t-5
        path = paths[i]
        plot!(ns, path[ns]; label="", lw=0.3, alpha=0.5)
    end
    for i in max(1, t-4):min(t, length(paths))
        path = paths[i]
        j = i - (t-5)
        plot!(ns, path[ns]; label="", lw=0.3+0.2j, alpha=min(1, 0.3+0.2j))
    end
    plot!(; kwargs...)
end

function gif_lln(; n=10^4, p=1/6, p_str="1/6", niters=100, ylim=(p-0.05, p+0.05), ytick=0:0.01:1, kwargs...)
    paths = [cumsum(rand(Bernoulli(p), n)) ./ (1:n) for _ in 1:niters]
    @time anim = @animate for t in [1:(niters+5); fill(niters+5, 10)]
        plot_paths(paths, t, p; xlim=(-0.02n, 1.02n), ylim)
        hline!([p]; label="p = $p_str", c=:black, alpha=0.7)
        plot!(range(0, n, 1001), x->p+2√(p*(1-p)/x); label="p±2√(p(1-p)/x)", c=:red)
        plot!(range(0, n, 1001), x->p-2√(p*(1-p)/x); label="", c=:red)
        plot!(; ytick, kwargs...)
        title!("n = $n,  p = $p_str,  t = $(min(niters, t))")
    end
    @time gif(anim, "lln$n.gif", fps=5)
end

Random.seed!(4649373)
gif_lln(; n=10^4, p=1/6, p_str="1/6", niters=100, ylim=(1/6-0.05, 1/6+0.05), ytick=0:0.01:1)

# %%
Random.seed!(4649373)
gif_lln(; n=10^5, p=1/6, p_str="1/6", niters=100, ylim=(1/6-0.016, 1/6+0.016), ytick=0:0.002:1)

# %%
Random.seed!(4649373)
gif_lln(; n=10^6, p=1/6, p_str="1/6", niters=100, ylim=(1/6-0.005, 1/6+0.005), ytick=0:0.001:1)

# %%
Random.seed!(4649373)
gif_lln(; n=10^7, p=1/6, p_str="1/6", niters=100, ylim=(1/6-0.0016, 1/6+0.0016), ytick=0:0.0002:1)

# %%
