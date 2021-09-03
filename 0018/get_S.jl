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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
using BenchmarkTools

# %%
W = trues(30000, 30001)
W[rand(1:length(W), round(Int, √length(W)))] .= 0
W

# %%
function get_S(W)
    zeroidxs = getindex.(findall(iszero, W), [1 2])
    S1 = Set(zeroidxs[:,1])
    S2 = Set(zeroidxs[:,2])
    return S1, S2
end

# %%
@btime get_S($W)

# %%
function get_S1(W)
    zeroidxs = getindex.(findall(iszero, W), [1 2])
    S1 = unique(zeroidxs[:,1])
    S2 = unique(zeroidxs[:,2])
    return S1, S2
end

# %%
@btime get_S1($W)

# %%
sort.(unique.(collect.(get_S(W)))) == sort.(get_S1(W))

# %%
function get_S2(W)
    m, n = size(W)
    b1, b2 = fill(false, m), fill(false, n)
    Threads.@threads for j in 1:n
        for i in 1:m
            @inbounds if W[i, j] == 0
                b1[i] = b2[j] = true
            end
        end
    end
    b1, b2
end

# %%
@btime get_S2($W)

# %%
sort.(get_S1(W)) == findall.(get_S2(W))

# %%
get_S(W) |> typeof

# %%
get_S1(W) |> typeof

# %%
get_S2(W) |> typeof

# %%
