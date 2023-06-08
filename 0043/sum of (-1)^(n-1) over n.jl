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
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# $\displaystyle\sum_{n=1}^L \frac{(-1)^{n-1}}{n}$

# %%
log(2)

# %%
a(n) = (-1)^(n-1)/n
L = 10^7
@time sum(a, 1:L)
@time sum(a, 1:L)
@time sum(a, 1:L)

# %%
using LoopVectorization

function f(L = 10^7)
    s = 0.0
    @turbo for n in 1:L
        s += (-1)^(n-1)/n
    end
    s
end

@time f(10^7)
@time f(10^7)
@time f(10^7)

# %%
