# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://twitter.com/dannchu/status/1443524638252810248

# %%
using Plots

function f(n, a=falses(6))
    s = 0
    for _ in 1:n
        i = rand(1:6)
        a[i] && continue
        a[i] = true
        s += i
        all(a) && break
    end
    s
end

for n in 1:10
    @time X = [f(n) for _ in 1:10^6]; flush(stdout)
    histogram(X; norm=true, alpha=0.3, bin=0.5:21.5, xtick=1:21, label="", title="n = $n") |> display
end

# %%
