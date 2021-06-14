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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %%
function pimc(n)
    c = 0
    for _ in 1:n
        c += rand()^2 + rand()^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc(10^9)
@time pimc(10^9)
@time pimc(10^9)

# %%
using Random: default_rng

function pimc2(n)
    rng = default_rng()
    c = 0
    for _ in 1:n
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc2(10^9)
@time pimc2(10^9)
@time pimc2(10^9)

# %%
using Random: MersenneTwister

function pimc3(n)
    rng = MersenneTwister()
    c = 0
    for _ in 1:n
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc3(10^9)
@time pimc3(10^9)
@time pimc3(10^9)

# %%
function pimc(n)
    c = 0
    for _ in 1:n
        c += rand()^2 + rand()^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc(10^9)
@time pimc(10^9)
@time pimc(10^9)

# %%
using Random: default_rng

function pimc2(n)
    rng = default_rng()
    c = 0
    for _ in 1:n
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc2(10^9)
@time pimc2(10^9)
@time pimc2(10^9)

# %%
using Random: MersenneTwister

function pimc3(n)
    rng = MersenneTwister()
    c = 0
    for _ in 1:n
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc3(10^9)
@time pimc3(10^9)
@time pimc3(10^9)

# %%
function pimc(n)
    c = 0
    for _ in 1:n
        c += rand()^2 + rand()^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc(10^9)
@time pimc(10^9)
@time pimc(10^9)

# %%
using Random: default_rng

function pimc2(n)
    rng = default_rng()
    c = 0
    for _ in 1:n
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc2(10^9)
@time pimc2(10^9)
@time pimc2(10^9)

# %%
using Random: MersenneTwister

function pimc3(n)
    rng = MersenneTwister()
    c = 0
    for _ in 1:n
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/n
end

@show VERSION
@time pimc3(10^9)
@time pimc3(10^9)
@time pimc3(10^9)

# %%
