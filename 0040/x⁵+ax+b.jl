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
#     display_name: Julia 1.9.0-beta2
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# * https://www.wolframalpha.com/input/?i=x%5E5%2Bax%2Bb%3D0
# * https://www.wolframalpha.com/input?i=x%5E5%2B3x%2B1%3D0

# %% [markdown]
# $x^5+ax+b=0$ の解の1つを
#
# $$x = -(b _4 F_3(1/5, 2/5, 3/5, 4/5 ;1/2, 3/4, 5/4 ;-(3125 b^4)/(256 a^5)))/a$$
#
# で作れる.

# %%
using HypergeometricFunctions

f(a, b, x) = x*(x^4 + a) + b

# x = -(b _4 F_3(1/5, 2/5, 3/5, 4/5 ;1/2, 3/4, 5/4 ;-(3125 b^4)/(256 a^5)))/a
x1(a, b) = -(b/a) * pFq([1//5, 2//5, 3//5, 4//5], [1//2, 3//4, 5//4], -(3125//256)*(b/a)^4/a)

function test_x1(a, b)
    @show a, b
    @show -(3125//256)*(b/a)^4/a
    @show x = x1(a, b)
    @show f(a, b, x)
end

# %%
test_x1(3, 1);

# %%
test_x1(big(3), big(1));

# %%
test_x1(1, 1);

# %%
test_x1(big(1), big(1));

# %%
test_x1(1, 3);

# %%
test_x1(big(1), big(3));

# %%
for i in 1:20
    z = test_x1(randn(BigFloat), randn(BigFloat))
    if abs(z) > 0.01
        println("誤差が大きい.")
    end
    println()
end

# %%
