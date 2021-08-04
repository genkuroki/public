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
mutable struct Mean_covar
    mu::Array{Float64,2}
    sigma::Array{Float64,2}
    w::Float64 
end

const B = Array{Mean_covar,1}(undef, 5)
display(B)

# %%
B[1] = Mean_covar(zeros(2,2), zeros(2,2), 0.0)
display(B)

# %%
