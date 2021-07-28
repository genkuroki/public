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
#     display_name: Julia 1.7.0-beta3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Statistics
using Plots
using Distributed

rmprocs(procs()[2:end])
@show addprocs()

function mcpi1(n)
    c = @distributed (+) for _ in 1:n
        rand()^2 + rand()^2 ≤ 1
    end
    4c/n
end

function mcpi2(n)
    c = @distributed (+) for _ in 1:n
        x, y = rand(), rand()
        if x + y < 1
            x, y = 1 - y, 1 - x
        end
        x^2 + y^2 ≤ 1
    end
    2c/n + 2
end

@time pi1 = [mcpi1(10^6) for _ in 1:10^4]
@time pi2 = [mcpi2(10^6) for _ in 1:10^4];

@show mean(pi1), std(pi1)
@show mean(pi2), std(pi2)

stephist(pi1; norm=true, label="mcpi1")
stephist!(pi2; norm=true, label="mcpi2")

# %%
