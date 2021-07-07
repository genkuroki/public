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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using Tullio, LoopVectorization

I, J, K = 400, 200, 100
f1 = rand(I, 1, K)
f2 = rand(I, J, 1)
dev = rand(I, 1, 1)

function make_R(f1, f2, dev)
    I, J, K = axes(dev, 1), axes(f2, 2), axes(f1, 3)
    R = similar(dev, I, J, K)
    for i in I, j in J, k in K
        @inbounds R[i, j, k] = f1[i, 1, k] * f2[i, j, 1] / dev[i, 1, 1]
    end
    R
end

R1 = make_R(f1, f2, dev)
@tullio R2[i, j, k] := f1[i, 1, k] * f2[i, j, 1] / dev[i, 1, 1]
R3 = @. f1 * f2 / dev
R4 = @turbo @. f1 * f2 / dev
R5 = @tturbo @. f1 * f2 / dev

@show R1 ≈ R2 ≈ R3 ≈ R4 ≈ R5
@show typeof(R1) size(R1);

# %%
using BenchmarkTools
@btime R1 = make_R($f1, $f2, $dev)
@btime @tullio R2[i, j, k] := $f1[i, 1, k] * $f2[i, j, 1] / $dev[i, 1, 1]
@btime R3 = @. $f1 * $f2 / $dev
@btime R4 = @turbo @. $f1 * $f2 / $dev
@btime R5 = @tturbo @. $f1 * $f2 / $dev
;

# %%
