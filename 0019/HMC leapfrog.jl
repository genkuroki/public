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
#     display_name: Julia 1.6.4
#     language: julia
#     name: julia-1.6
# ---

# %%
# https://twitter.com/iitenki_moruten/status/1434070751884177415

using Plots
using StatsBase
using Random
using Printf
using Distributions
using BenchmarkTools
using QuadGK
using SpecialFunctions
using Zygote

# Generating Gaussian distribution (1-dim) by using HMC method
# 改良後

v(x) = 0.5*x^2
u(x) = x^2*(x-3)^2/2

function S(x,n_size)
    sum(u(x[i]) for i=1:n_size)
end

function dSdx(x,k)
    return u'(x) 
end

function H(x,p,n_size)
    K = sum(p[i]^2 for i=1:n_size)
    return S(x,n_size) + 0.5*K
end

#１変数更新のHMC法
function HMC_1dim(x,parameter)
        n_size,n_tau,τ = parameter
        #初期化
        p_ini = rand(Normal(0,1),n_size)
        p_fin = zeros(n_size)
        x_ini = zeros(n_size)
        x_fin = zeros(n_size)
    
        #初期配位を保存
        for i=1:n_size x_ini[i] = x[i] end
    
        #初期配位のハミルトニアンを計算
        H_ini = H(x_ini,p_ini,n_size)
    
        #リープフロッグで時間発展
        leapfrog!(x,p_ini,parameter)
    
        #時間発展後の配位を保存
        for i=1:n_size x_fin[i] = x[i]; p_fin[i] = p_ini[i] end
    
        #時間発展後のハミルトニアンを計算
        H_fin = H(x_fin,p_fin,n_size)
    
        #Metropolis check
        r=rand()
        if r<exp(H_ini - H_fin)
            for i=1:n_size x[i] = x_fin[i] end
        else
            for i=1:n_size x[i] = x_ini[i] end
        end
end

#x,pをリープフロッグ法で更新する．
function leapfrog!(x,p,parameter)
    n_size,n_tau,τ = parameter
    Δτ  = τ/n_tau 
    p0  = 0.0
    p1  = 0.0
    x05 = 0.0
    x15 = 0.0
    for k=1:n_size
        #initial setting : τ0
        p0  = p[k]
        x05 = x[k] + p0 * 0.5*Δτ
    
        #(n-1)step : τ0 -> τn-1 
        for i=1:n_tau-1
            p1  = p0 - dSdx(x05,k)*Δτ
            x15 = x05 + p1*Δτ
            p0  = p1
            x05 = x15
        end
    
        #final step : τn-1 -> τn
        p_end = p1 - dSdx(x15,k)*Δτ
        x_end = x15 + p_end*0.5*Δτ
        
        x[k]  = x_end
        p[k]  = p_end
    end
end

function gen_configs(Ω,n_con,parameter)
    n_size,n_tau,τ = parameter
    x = zeros(n_size)
    for j=1:n_con
        HMC_1dim(x,parameter)
        Ω[j] = x[1]
    end
end

##########################################################

"""q = [n_size,n_tau,τ]"""

#パラメータベクトル
n_config = 1_000_000          # # of configurations
n_size   = 1
n_tau    = 100
τ        = 5
q 　　　  = [n_size,n_tau,τ]
Ω        = zeros(n_config)

#result = @benchmark gen_configs(n_config,q)
@time gen_configs(Ω,n_config,q)

##########################################################


###plot###
gauss(x) = exp(-0.5*x^2)/sqrt(2*pi)
pdf_u(x) = exp(-u(x))
Z0,error = quadgk(pdf_u,-Inf,Inf)
normed_u(x) = pdf_u(x)/Z0

xr = range(-3,6,length=1000)

plot(Ω,st=:histogram,nbins=50,norm=:pdf,alpha=0.3,label="HMC's sample",legend=:topleft)
plot!(xr,normed_u.(xr),label="pdf_ana")

# %%
##########################################################

"""q = [n_size,n_tau,τ]"""

#パラメータベクトル
n_config = 1_000_000          # # of configurations
n_size   = 1
n_tau    = 100
τ        = 5
q 　　　  = [n_size,n_tau,τ]
Ω        = zeros(n_config)

#result = @benchmark gen_configs(n_config,q)
@time gen_configs(Ω,n_config,q)

##########################################################


###plot###
gauss(x) = exp(-0.5*x^2)/sqrt(2*pi)
pdf_u(x) = exp(-u(x))
Z0,error = quadgk(pdf_u,-Inf,Inf)
normed_u(x) = pdf_u(x)/Z0

xr = range(-3,6,length=1000)

plot(Ω,st=:histogram,nbins=50,norm=:pdf,alpha=0.3,label="HMC's sample",legend=:topleft)
plot!(xr,normed_u.(xr),label="pdf_ana")

# %% [markdown]
# 以下のコードでの主な変更点
#
# * `r < exp(H_ini - H_fin)` → `r ≤ min(1, exp(H_ini - H_fin))` (Metropolis-Hastings化)
# * `n_size` を函数内から削除 (必要なら `n_size = length(x)` とできる)
# * `parameter` に作業用配列 `p_ini`, `p_fin`, `x_ini`, `x_fin` を含めて、メモリアロケーションを削減。
# * `parameter` にポテンシャル函数 `u` も含め、函数 `H` などに渡すようにした。
# * Zygote.jl ではなく、ForwardDiff.jl を使うようにして、初回実行時の遅延の短縮。
#
# References
#
# * https://docs.julialang.org/en/v1/manual/performance-tips/#Pre-allocating-outputs
# * https://docs.julialang.org/en/v1/manual/performance-tips/#More-dots:-Fuse-vectorized-operations
# * https://github.com/genkuroki/public/blob/main/0018/HMC%20leapfrog.ipynb

# %%
using Plots
# using StatsBase
# using Random
# using Printf
using Distributions
using BenchmarkTools
using QuadGK
# using SpecialFunctions

module Revised

using Random
using Distributions
using ForwardDiff

# Generating Gaussian distribution (1-dim) by using HMC method
# 改良後

S(x, parameter) = ((u,) = parameter; sum(u, x))
dSdx(x, k, parameter) = ((u,) = parameter; ForwardDiff.derivative(u, x))
square(x) = x^2
H(x, p, parameter) = S(x, parameter) + 0.5sum(square, p)

#１変数更新のHMC法
function HMC_1dim!(x, parameter)
        _, n_tau, τ, p_ini, p_fin, x_ini, x_fin = parameter
        
        #初期化
        rand!(Normal(0,1), p_ini)
    
        #初期配位を保存
        x_ini .= x
    
        #初期配位のハミルトニアンを計算
        H_ini = H(x_ini, p_ini, parameter)
    
        #リープフロッグで時間発展
        leapfrog!(x, p_ini, parameter)
    
        #時間発展後の配位を保存
        x_fin .= x
        p_fin .= p_ini
    
        #時間発展後のハミルトニアンを計算
        H_fin = H(x_fin, p_fin, parameter)
    
        #Metropolis check
        r = rand()
        if r ≤ min(1, exp(H_ini - H_fin))
            x .= x_fin
        else
            x .= x_ini
        end
end

#x,pをリープフロッグ法で更新する．
function leapfrog!(x, p, parameter)
    n_size = length(x)
    _, n_tau, τ = parameter
    Δτ  = τ/n_tau
    for k=1:n_size
        #initial setting : τ0
        p0  = p[k]
        x05 = x[k] + p0 * 0.5*Δτ
    
        #(n-1)step : τ0 -> τn-1
        local p1, x15
        for i=1:n_tau-1
            p1  = p0 - dSdx(x05, k, parameter)*Δτ
            x15 = x05 + p1*Δτ
            p0  = p1
            x05 = x15
        end
    
        #final step : τn-1 -> τn
        p_end = p1 - dSdx(x15, k, parameter)*Δτ
        x_end = x15 + p_end*0.5*Δτ
        
        x[k]  = x_end
        p[k]  = p_end
    end
end

function gen_configs(Ω, x, parameter)
    n_tau, τ = parameter
    for j in eachindex(Ω)
        HMC_1dim!(x, parameter)
        Ω[j] = x[1]
    end
end

end

v(x) = 0.5x^2
u(x) = x^2*(x-3)^2/2

##########################################################

n_config = 10^6          # # of configurations
n_size   = 1
n_tau    = 100
τ        = 5

Ω        = zeros(n_config)
x        = zeros(n_size)
p_ini    = similar(x)
p_fin    = similar(x)
x_ini    = similar(x)
x_fin    = similar(x)

"""q = (u, n_tau, τ, p_ini, p_fin, x_ini, x_fin)"""
q 　　　  = (u, n_tau, τ, p_ini, p_fin, x_ini, x_fin)

#result = @benchmark gen_configs(n_config,q)
@time Revised.gen_configs(Ω, x, q)

##########################################################

###plot###
gauss(x) = exp(-0.5*x^2)/sqrt(2*pi)
pdf_u(x) = exp(-u(x))
Z0,error = quadgk(pdf_u,-Inf,Inf)
normed_u(x) = pdf_u(x)/Z0

xr = range(-3,6,length=1000)

plot(Ω,st=:histogram,nbins=50,norm=:pdf,alpha=0.3,label="HMC's sample",legend=:topleft)
plot!(xr,normed_u.(xr),label="pdf_ana")

# %%
##########################################################

n_config = 10^6          # # of configurations
n_size   = 1
n_tau    = 100
τ        = 5

Ω        = zeros(n_config)
x        = zeros(n_size)
p_ini    = similar(x)
p_fin    = similar(x)
x_ini    = similar(x)
x_fin    = similar(x)

"""q = (u, n_tau, τ, p_ini, p_fin, x_ini, x_fin)"""
q 　　　  = (u, n_tau, τ, p_ini, p_fin, x_ini, x_fin)

#result = @benchmark gen_configs(n_config,q)
@time Revised.gen_configs(Ω, x, q)

##########################################################

###plot###
gauss(x) = exp(-0.5*x^2)/sqrt(2*pi)
pdf_u(x) = exp(-u(x))
Z0,error = quadgk(pdf_u,-Inf,Inf)
normed_u(x) = pdf_u(x)/Z0

xr = range(-3,6,length=1000)

plot(Ω,st=:histogram,nbins=50,norm=:pdf,alpha=0.3,label="HMC's sample",legend=:topleft)
plot!(xr,normed_u.(xr),label="pdf_ana")

# %%
