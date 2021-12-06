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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# https://twitter.com/iitenki_moruten/status/1467683477474930688

# %%
# Original: https://github.com/moruten/julia-code/blob/1dafca3e1a4e3b36445c2263e440f6e4056b90aa/2021-12-6-test-no1.ipynb

using Plots
using Random
using Distributions
using QuadGK

#============関数定義============================================#
function action(x)
    0.5*x*x
end

function deriv_action(x)
    x
end

function hamiltonian(x,p)
    action(x) + 0.5*p*p
end

function HMC_update(x,Nt,dt)
    #backup
    x_old = x
    p     = rand(Normal(0,1))
    
    #check
    H_ini = hamiltonian(x,p)
    x     = molecular_dynamics!(x,p,Nt,dt)
    H_fin = hamiltonian(x,p)
    
    r  = rand()
    ΔH = H_fin-H_ini
    if r < exp(-ΔH)
        return x,1 #accept
    else 
        return x_old,0
    end  
end

function molecular_dynamics!(x,p,Nt,dt) 
    force = 0.0
    #1/2step
        x = x + p * 0.5*dt
    #1~Nt-1 step
        for j=1:(Nt-1)
            force = deriv_action(x)
            p = p - dt * force
            x = x + dt * p
        end
    #Nt step
        p = p - dt * deriv_action(x)
        x = x + 0.5 * dt * p
        return x
end
#============関数終わり============================================#

#============計算=========================================================================#
  #セットアップ
    Ntest        = 300000
    Nt           = 20
    dt           = 1.0/Nt
    conf_vec     = zeros(Ntest)
    accept_count = 0
    ret          = 0
    x            = 0.0

    sumxx        = 0.0
    sumx         = 0.0

  #計算
    for i=1:Ntest
        x,ret         = HMC_update(x,Nt,dt)
        accept_count += ret
        conf_vec[i]   = x
    
        sumx         += x
        sumxx        += x*x
    end

    println("P(accept) = $(accept_count/Ntest)")
    println("<x>       = $(sumx/Ntest)")
    println("<xx>      = $(sumxx/Ntest)")

#=======確認=============================================================================#
    xr = range(-5,5,length=1000)
    f1(x) = exp(-0.5*x^2)
    f2(x) = exp(-x^2)
    Z1,error1 = quadgk(f1,-Inf,Inf)
    Z2,error2 = quadgk(f2,-Inf,Inf)
    g1(x) = f1(x)/Z1
    g2(x) = f2(x)/Z2
    histogram(conf_vec,norm=:true,label="data")
    plot!(xr,[g1.(xr) g2.(xr)],lw=3,label=["exp(-0.5x^2)/Z1" "exp(-x^2)/Z2"])
#========================================================================================#

# %%
# Revised 1

using Plots
using Random
using Distributions
using QuadGK

#============関数定義============================================#
function action(x)
    0.5*x*x
end

function deriv_action(x)
    x
end

function hamiltonian(x,p)
    action(x) + 0.5*p*p
end

function HMC_update(x,Nt,dt)
    #backup
    x_old = x
    p     = rand(Normal(0,1))
    
    #check
    H_ini = hamiltonian(x,p)
    x, p = molecular_dynamics!(x,p,Nt,dt)    # <========== ここ
    H_fin = hamiltonian(x,p)
    
    r  = rand()
    ΔH = H_fin-H_ini
    if r < exp(-ΔH)
        return x,1 #accept
    else 
        return x_old,0
    end  
end

function molecular_dynamics!(x,p,Nt,dt) 
    force = 0.0
    #1/2step
        x = x + p * 0.5*dt
    #1~Nt-1 step
        for j=1:(Nt-1)
            force = deriv_action(x)
            p = p - dt * force
            x = x + dt * p
        end
    #Nt step
        p = p - dt * deriv_action(x)
        x = x + 0.5 * dt * p
        return x, p                          # <========== ここ
end
#============関数終わり============================================#

#============計算=========================================================================#
  #セットアップ
    Ntest        = 300000
    Nt           = 20
    dt           = 1.0/Nt
    conf_vec     = zeros(Ntest)
    accept_count = 0
    ret          = 0
    x            = 0.0

    sumxx        = 0.0
    sumx         = 0.0

  #計算
    for i=1:Ntest
        x,ret         = HMC_update(x,Nt,dt)
        accept_count += ret
        conf_vec[i]   = x
    
        sumx         += x
        sumxx        += x*x
    end

    println("P(accept) = $(accept_count/Ntest)")
    println("<x>       = $(sumx/Ntest)")
    println("<xx>      = $(sumxx/Ntest)")

#=======確認=============================================================================#
    xr = range(-5,5,length=1000)
    f1(x) = exp(-0.5*x^2)
    f2(x) = exp(-x^2)
    Z1,error1 = quadgk(f1,-Inf,Inf)
    Z2,error2 = quadgk(f2,-Inf,Inf)
    g1(x) = f1(x)/Z1
    g2(x) = f2(x)/Z2
    histogram(conf_vec,norm=:true,label="data")
    plot!(xr,[g1.(xr) g2.(xr)],lw=3,label=["exp(-0.5x^2)/Z1" "exp(-x^2)/Z2"])
#========================================================================================#

# %%
# Revised 2

using Plots
using Random
using Distributions
using QuadGK

#============関数定義============================================#
function action(x)
    x^2/2
end

function deriv_action(x)
    x
end

function hamiltonian(x,p)
    action(x) + 0.5*p*p
end

function HMC_update(x,Nt,dt)
    #backup
    x_old = x
    p     = rand(Normal(0,1))
    
    #check
    H_ini = hamiltonian(x,p)
    x, p = molecular_dynamics!(x,p,Nt,dt)     # <========== ここ
    H_fin = hamiltonian(x,p)
    
    r  = rand()
    ΔH = H_fin-H_ini
    if r < exp(-ΔH)
        return x,1 #accept
    else 
        return x_old,0
    end  
end

function molecular_dynamics!(x,p,Nt,dt)
    p -= deriv_action(x) * dt/2
    x += p * dt
    for j in 2:Nt
        p -= deriv_action(x) * dt
        x += p * dt
    end
    p -= deriv_action(x) * dt/2
    return x, p                    # <========== ここ
end
#============関数終わり============================================#

#============計算=========================================================================#
  #セットアップ
    Ntest        = 300000
    Nt           = 20
    dt           = 1.0/Nt
    conf_vec     = zeros(Ntest)
    accept_count = 0
    ret          = 0
    x            = 0.0

    sumxx        = 0.0
    sumx         = 0.0

  #計算
    for i=1:Ntest
        x,ret         = HMC_update(x,Nt,dt)
        accept_count += ret
        conf_vec[i]   = x
    
        sumx         += x
        sumxx        += x*x
    end

    println("P(accept) = $(accept_count/Ntest)")
    println("<x>       = $(sumx/Ntest)")
    println("<xx>      = $(sumxx/Ntest)")

#=======確認=============================================================================#
    xr = range(-2,2,length=1000)
    f(x) = exp(-action(x))
    Z,error1 = quadgk(f,-Inf,Inf)
    g(x) = f(x)/Z
    histogram(conf_vec,norm=:true,label="data")
    plot!(xr,g1.(xr),lw=3,label="exp(-action(x))/Z", legend=:outertop)
#========================================================================================#

# %%
# Revised 1 - second test

using Plots
using Random
using Distributions
using QuadGK

#============関数定義============================================#
function action(x)
    3(x^2 - 1)^2
end

function deriv_action(x)
    6x*(x^2 - 1)
end

function hamiltonian(x,p)
    action(x) + 0.5*p*p
end

function HMC_update(x,Nt,dt)
    #backup
    x_old = x
    p     = rand(Normal(0,1))
    
    #check
    H_ini = hamiltonian(x,p)
    x, p = molecular_dynamics!(x,p,Nt,dt)     # <========== ここ
    H_fin = hamiltonian(x,p)
    
    r  = rand()
    ΔH = H_fin-H_ini
    if r < exp(-ΔH)
        return x,1 #accept
    else 
        return x_old,0
    end  
end

function molecular_dynamics!(x,p,Nt,dt) 
    force = 0.0
    #1/2step
        x = x + p * 0.5*dt
    #1~Nt-1 step
        for j=1:(Nt-1)
            force = deriv_action(x)
            p = p - dt * force
            x = x + dt * p
        end
    #Nt step
        p = p - dt * deriv_action(x)
        x = x + 0.5 * dt * p
        return x, p                    # <========== ここ
end
#============関数終わり============================================#

#============計算=========================================================================#
  #セットアップ
    Ntest        = 300000
    Nt           = 20
    dt           = 1.0/Nt
    conf_vec     = zeros(Ntest)
    accept_count = 0
    ret          = 0
    x            = 0.0

    sumxx        = 0.0
    sumx         = 0.0

  #計算
    for i=1:Ntest
        x,ret         = HMC_update(x,Nt,dt)
        accept_count += ret
        conf_vec[i]   = x
    
        sumx         += x
        sumxx        += x*x
    end

    println("P(accept) = $(accept_count/Ntest)")
    println("<x>       = $(sumx/Ntest)")
    println("<xx>      = $(sumxx/Ntest)")

#=======確認=============================================================================#
    xr = range(-2,2,length=1000)
    f(x) = exp(-action(x))
    Z,error1 = quadgk(f,-Inf,Inf)
    g(x) = f(x)/Z
    histogram(conf_vec,norm=:true,label="data")
    plot!(xr,g.(xr),lw=3,label="exp(-action(x))/Z", legend=:outertop)
#========================================================================================#

# %%
