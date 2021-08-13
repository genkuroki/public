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
"""https://github.com/moruten/julia-code/blob/99b50bf329195ba102a9496602a22138f07d186b/Test-2021-8-12-no6.ipynb"""
module Original

using Plots
using Random

function any_width_rand!(n_rand,width_l,width_r)
    rands=zeros(n_rand)
    L=width_r - width_l
    for i=1:n_rand
        rands[i] = width_l + L*rand()
    end
    return rands
end

function action_cal(x)
    exp(0.5*x^2)
end

function MP_test(x_ini,x_pro)
    S_ini = action_cal(x_ini)
    S_pro = action_cal(x_pro)
    
    if exp(S_ini - S_pro) > rand()
        return x_pro
    else 
        return x_ini
    end
end

n_rand=1_000_000

MP_rands=zeros(n_rand)

x_propose  = 0.0
x_ini      = 0.0
width_rand = 2.0

@time for i=2:n_rand
    x_propose   = x_ini + any_width_rand!(1,-width_rand,width_rand)[1]
    MP_rands[i] = MP_test(x_ini,x_propose)
    global x_ini       = MP_rands[i]
end

f(x) = exp(-0.5*x^2)/sqrt(2*pi)
xr   = range(-3,3,length=10001)

plot(MP_rands,st=:histogram,norm=:true,nbins=100,label="Metropolis method")
plot!(xr,f.(xr),lw=3,label="Gaussian") |> display

end

# %%
"""https://github.com/moruten/julia-code/blob/99b50bf329195ba102a9496602a22138f07d186b/Test-2021-8-12-no6.ipynb"""
module Original_plus

using Plots
using Random

function any_width_rand!(n_rand,width_l,width_r)
    rands=zeros(n_rand)
    L=width_r - width_l
    for i=1:n_rand
        rands[i] = width_l + L*rand()
    end
    return rands
end

function action_cal(x)
    exp(0.5*x^2)
end

function MP_test(x_ini,x_pro)
    S_ini = action_cal(x_ini)
    S_pro = action_cal(x_pro)
    
    if exp(S_ini - S_pro) > rand()
        return x_pro
    else 
        return x_ini
    end
end

n_rand=1_000_000

MP_rands=zeros(n_rand)

x_propose  = 0.0
x_ini      = 0.0
width_rand = 2.0

@time for i=2:n_rand
    x_propose   = x_ini + any_width_rand!(1,-width_rand,width_rand)[1]
    MP_rands[i] = MP_test(x_ini,x_propose)
    global x_ini       = MP_rands[i]
end

f(x) = exp(-0.5*x^2)/sqrt(2*pi)
xr   = range(-3,3,length=10001)

using QuadGK
G(x) = exp(-exp(0.5*x^2))
Z = quadgk(G, -Inf, Inf)[1]
g(x) = G(x)/Z

histogram(MP_rands, norm=true, nbins=100, label="Metropolis method", alpha=0.3)
plot!(xr, f.(xr), lw=3, label="Gaussian")
plot!(xr, g.(xr), lw=3, label="exp(-exp(0.5*x^2))/Z", ls=:dash) |> display

end

# %%
module Corrected

using Plots
using Random

function any_width_rand!(n_rand,width_l,width_r)
    rands=zeros(n_rand)
    L=width_r - width_l
    for i=1:n_rand
        rands[i] = width_l + L*rand()
    end
    return rands
end

function action_cal(x)
    exp(0.5*x^2)
end

function MP_test(x_ini,x_pro)
    S_ini = action_cal(x_ini)
    S_pro = action_cal(x_pro)
    
    if S_ini / S_pro > rand() # <-- 変更はこの行だけ！
        return x_pro
    else 
        return x_ini
    end
end

n_rand=1_000_000

MP_rands=zeros(n_rand)

x_propose  = 0.0
x_ini      = 0.0
width_rand = 2.0

@time for i=2:n_rand
    x_propose   = x_ini + any_width_rand!(1,-width_rand,width_rand)[1]
    MP_rands[i] = MP_test(x_ini,x_propose)
    global x_ini       = MP_rands[i]
end

f(x) = exp(-0.5*x^2)/sqrt(2*pi)
xr   = range(-5,5,length=10001)

plot(MP_rands,st=:histogram,norm=:true,nbins=100,label="Metropolis method", alpha=0.3)
plot!(xr,f.(xr),lw=3,label="Gaussian") |> display

end

# %%
module MH

using Plots
using Random

function any_width_rand!(n_rand,width_l,width_r)
    rands=zeros(n_rand)
    L=width_r - width_l
    for i=1:n_rand
        rands[i] = width_l + L*rand()
    end
    return rands
end

function action_cal(x)
    exp(0.5*x^2)
end

function MP_test(x_ini,x_pro)
    S_ini = action_cal(x_ini)
    S_pro = action_cal(x_pro)
    a = S_ini/S_pro
    if a > 1 || rand() ≤ a
        return x_pro
    else 
        return x_ini
    end
end

n_rand=1_000_000

MP_rands=zeros(n_rand)

x_propose  = 0.0
x_ini      = 0.0
width_rand = 2.0

@time for i=2:n_rand
    x_propose   = x_ini + any_width_rand!(1,-width_rand,width_rand)[1]
    MP_rands[i] = MP_test(x_ini,x_propose)
    global x_ini       = MP_rands[i]
end

f(x) = exp(-0.5*x^2)/sqrt(2*pi)
xr   = range(-5,5,length=10001)

plot(MP_rands,st=:histogram,norm=:true,nbins=100,label="MH method", alpha=0.3)
plot!(xr,f.(xr),lw=3,label="Gaussian") |> display

end

# %%
using Plots
using Distributions

function update(H, prop, x)
    y = prop(x)
    a = exp(H(x) - H(y))
    a > 1 || rand() ≤ a ? y : x
end

function iter_update(H, prop, x0, N; warmup=min(N ÷ 10, 10^3))
    x = x0
    for _ in 1:warmup
        x = update(H, prop, x)
    end
    X = similar([x], N)
    X[1] = x
    for i in 2:N
        X[i] = update(H, prop, X[i-1])
    end
    X
end

H(x) = x^2/2
prop(x) = x + 0.2randn()
x0 = randn()
N = 10^6
X = iter_update(H, prop, x0, N)

histogram(X; norm=true, alpha=0.3, bin=100, label="MH sample")
plot!(x -> pdf(Normal(), x), -5, 5; ls=:dash, lw=2, label="normal dist.")

# %%
H(x) = x ≤ 0 ? Inf : x
prop(x) = x + 0.1randn()
x0 = abs(randn())
N = 10^6
X = iter_update(H, prop, x0, N)

histogram(X; norm=true, alpha=0.3, bin=100, label="MH sample")
plot!(x -> pdf(Exponential(), x), eps(), 10; label="exp. dist.")
plot!(xlim=(-0.1, 8))

# %%
H(x) = x ≤ 0 ? Inf : x/3 - (4 - 1)*log(x)
prop(x) = x + randn()
x0 = abs(randn())
N = 10^6
X = iter_update(H, prop, x0, N)

histogram(X; norm=true, alpha=0.3, bin=100, label="MH sample")
plot!(x -> pdf(Gamma(4, 3), x), eps(), 50; label="Gamma(4, 3)")

# %%
