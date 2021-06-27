# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.7.0-beta2
#     language: julia
#     name: julia-1.7
# ---

# %%
function f(n)
    s = 0.0
    for (i1, i2, i3, i4, i5) in Iterators.product(0:n-1, 0:n-1, 0:n-1, 0:n-1, 0:n-1)
        x = evalpoly(n, (i5, i4, i3, i2, i1)) + 1
        s += sin(x)/x
    end
    2s + 1
end

function g(n)
    s = 0.0
    for i1 in 0:n-1, i2 in 0:n-1, i3 in 0:n-1, i4 in 0:n-1, i5 in 0:n-1
        x = evalpoly(n, (i1, i2, i3, i4, i5)) + 1
        s += sin(x)/x
    end
    2s + 1
end

function h(n)
    s = 0.0
    for i in 0:n^5-1
        x = float(i) + 1
        s += sin(x)/x
    end
    2s + 1
end

using BenchmarkTools
@show VERSION
@show f(20)
@show g(20)
@show h(20)
@show f(20) ≈ g(20) ≈ h(20)
@btime f(20)
@btime g(20)
@btime h(20);

# %%
function f(n)
    s = 0.0
    for (i1, i2, i3, i4, i5) in Iterators.product(0:n-1, 0:n-1, 0:n-1, 0:n-1, 0:n-1)
        x = evalpoly(n, (i5, i4, i3, i2, i1)) + 1
        s += sin(x)/x
    end
    2s + 1
end

function g(n)
    s = 0.0
    for i1 in 0:n-1, i2 in 0:n-1, i3 in 0:n-1, i4 in 0:n-1, i5 in 0:n-1
        x = evalpoly(n, (i1, i2, i3, i4, i5)) + 1
        s += sin(x)/x
    end
    2s + 1
end

function h(n)
    s = 0.0
    for i in 0:n^5-1
        x = float(i) + 1
        s += sin(x)/x
    end
    2s + 1
end

using BenchmarkTools
@show VERSION
@show f(20)
@show g(20)
@show h(20)
@show f(20) ≈ g(20) ≈ h(20)
@btime f(20)
@btime g(20)
@btime h(20);

# %%
