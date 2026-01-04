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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# https://colab.research.google.com/github/genkuroki/public/blob/main/0055/happy%20numbers.ipynb

# %%
struct IsHappy{T<:Integer} ishappy::Bool; seq::Vector{T} end

function Base.show(io::IO, x::IsHappy)
    (; ishappy, seq) = x
    n = first(seq)
    print(io, n, " is ")
    ishappy ? print(io, "happy: ") : print(io, "unhappy: ")
    print(io, n)
    for k in @view(seq[begin+1:end]) print(io, '→', k) end
end

function IsHappy(n::T) where T<:Integer
    @assert n > 0
    seq = [n]
    while true
        n == 1 && break
        n = sum(k -> k^2, digits(T, n))
        n ∈ seq && (push!(seq, n); break)
        push!(seq, n)
    end
    ishappy = n == 1
    IsHappy(ishappy, seq)
end

IsHappy(2026)

# %%
IsHappy(2026) |> dump

# %%
IsHappy(Int32(2026)) |> dump

# %%
IsHappy(big(2026)) |> dump

# %%
IsHappy(1999)

# %%
ENV["LINES"] = 1000
[IsHappy(n) for n in 1900:2100]

# %%
@code_warntype IsHappy(1234)

# %%
