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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
powermod(10, 10, 2020)

# %%
let y = 2020, b = 0
    for x in 0:98
        b += (powermod(10, 99, y) + powermod(10, x, y)) % y == 0
    end
    b
end

# %%
let y = 2020
    sum((powermod(10, 99, y) + powermod(10, x, y)) % y == 0 for x in 0:98)
end

# %%
let y = 2020
    sum(x -> (powermod(10, 99, y) + powermod(10, x, y)) % y == 0, 0:98)
end

# %%
let y = 2020
    count((powermod(10, 99, y) + powermod(10, x, y)) % y == 0 for x in 0:98)
end

# %%
let y = 2020
    count(x -> (powermod(10, 99, y) + powermod(10, x, y)) % y == 0, 0:98)
end

# %%
