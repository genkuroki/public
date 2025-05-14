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
#     display_name: Julia current stable release
#     language: julia
#     name: julia
# ---

# %% [markdown]
# https://x.com/dannchu/status/1922545316819014000

# %%
versioninfo()

# %%
using BenchmarkTools

function row_major()
    n = 20_000
    x = randn(n, n)
    y = similar(x)
    @inbounds for r in 1:n
        for c in 1:n
            y[r, c] = max(0, x[r, c])
        end
    end
    return nothing
end

function col_major()
    n = 20_000
    x = randn(n, n)
    y = similar(x)
    @inbounds for c in 1:n
        for r in 1:n
            y[r, c] = max(0, x[r, c])
        end
    end
    return nothing
end

function use_eachindex()
    n = 20_000
    x = randn(n, n)
    y = similar(x)
    @inbounds for i in eachindex(x)
        y[i] = max(0, x[i])
    end
    return nothing
end

function max0dot()
    n = 20_000
    x = randn(n, n)
    y = max.(0, x)
    return nothing
end

println("Benchmarking...")

@btime row_major()
@btime col_major()
@btime use_eachindex()
@btime max0dot()

# %%
