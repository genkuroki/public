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

# %% [markdown]
# https://discourse.julialang.org/t/concrete-example-0805-how-to-prevent-using-global-variable-in-the-module/65917

# %%
module Mixture

export musigma

mutable struct Mean_covar
    mu::Array{Float64,2}
    sigma::Array{Float64,2}
    w::Float64
end

global const musigma = Array{Mean_covar,1}()  

function mean_covar_init(kmix::Int64,dim_p::Int64
                        ,weight::Array{Float64,1}
                        ,sigma::Array{Float64,3}
                        ,mu::Array{Float64,2})
    @assert length(weight) == kmix
    @assert size(sigma) == (kmix,dim_p,dim_p)
    @assert size(mu) == (kmix,dim_p)
    resize!(musigma, kmix) 
    for k in 1:kmix
        musigma[k] = Mean_covar(zeros(dim_p,1),zeros(dim_p,dim_p),0.0)
        musigma[k].mu[1,1] = mu[k,1]
        musigma[k].mu[2,1] = mu[k,2]
        musigma[k].sigma[1,1] = sigma[k,1,1]
        musigma[k].sigma[2,2] = sigma[k,2,2]
        musigma[k].w = weight[k]
    end
    return nothing
end

end

kmix = 5
dim_p = 3
weight = rand(kmix)
sigma = rand(kmix, dim_p, dim_p)
mu = rand(kmix, dim_p)

using .Mixture
Mixture.mean_covar_init(kmix, dim_p, weight, sigma, mu)
@show typeof(musigma) size(musigma)
@show musigma[1].mu musigma[1].sigma musigma[1].w;

# %%
module O

export musigma

struct Mean_covar{T}
    mu::Vector{T}
    sigma::Matrix{T}
    w::Array{T,0}
end

function mean_covar_init(kmix, dim_p, 
        weight::AbstractVector{T},
        sigma::AbstractArray{T,3},
        mu::AbstractMatrix{T}) where T
    @assert size(weight) == (kmix,)
    @assert size(sigma) == (kmix, dim_p, dim_p)
    @assert size(mu) == (kmix, dim_p)
    musigma = Vector{Mean_covar{T}}(undef, kmix)
    for k in 1:kmix
        musigma[k] = Mean_covar(zeros(T, dim_p), zeros(T, dim_p, dim_p), fill(zero(T)))
        musigma[k].mu[1,1] = mu[k,1]
        musigma[k].mu[2,1] = mu[k,2]
        musigma[k].sigma[1,1] = sigma[k,1,1]
        musigma[k].sigma[2,2] = sigma[k,2,2]
        musigma[k].w[] = weight[k]
    end
    musigma
end

end

musigma64 = O.mean_covar_init(kmix, dim_p, weight, sigma, mu)
@show typeof(musigma64) size(musigma64)
@show musigma64[1].mu musigma64[1].sigma musigma64[1].w;

# %%
weight32 = rand(Float32, kmix)
sigma32 = rand(Float32, kmix, dim_p, dim_p)
mu32 = rand(Float32, kmix, dim_p)

musigma32 = O.mean_covar_init(kmix, dim_p, weight32, sigma32, mu32)
@show typeof(musigma32) size(musigma32)
@show musigma32[1].mu musigma32[1].sigma musigma32[1].w;

# %%
module P

export musigma

struct Mean_covar{M,S,W}
    mu::M
    sigma::S
    w::W
end

function mean_covar_init(kmix, dim_p, 
        weight::AbstractVector{T},
        sigma::AbstractArray{T,3},
        mu::AbstractMatrix{T}) where T
    @assert size(weight) == (kmix,)
    @assert size(sigma) == (kmix, dim_p, dim_p)
    @assert size(mu) == (kmix, dim_p)
    M = typeof(similar(mu, (0,)))
    S = typeof(similar(sigma, (0, 0)))
    W = typeof(similar(weight, ()))
    musigma = Vector{Mean_covar{M,S,W}}(undef, kmix)
    for k in 1:kmix
        musigma[k] = Mean_covar(zeros(T, dim_p), zeros(T, dim_p, dim_p), fill(zero(T)))
        musigma[k].mu[1,1] = mu[k,1]
        musigma[k].mu[2,1] = mu[k,2]
        musigma[k].sigma[1,1] = sigma[k,1,1]
        musigma[k].sigma[2,2] = sigma[k,2,2]
        musigma[k].w[] = weight[k]
    end
    musigma
end

end

musigma64 = P.mean_covar_init(kmix, dim_p, weight, sigma, mu)
@show typeof(musigma64) size(musigma64)
@show musigma64[1].mu musigma64[1].sigma musigma64[1].w;

# %%
musigma32 = O.mean_covar_init(kmix, dim_p, weight32, sigma32, mu32)
@show typeof(musigma32) size(musigma32)
@show musigma32[1].mu musigma32[1].sigma musigma32[1].w;

# %%
