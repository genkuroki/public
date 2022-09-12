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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %%
using BenchmarkTools
using SpecialFunctions

incbeta_I(x, a, b) = beta_inc(a, b, x)[1]

function ibetar_st(x::T, a, b;
        HUGE_VAL = T(Inf), FDeps = eps(T), rtol = eps(T), MAXITERS=200
    ) where T<:AbstractFloat
    a ≤ 0 && return HUGE_VAL
    if b ≤ 0
        x < 1 && return zero(T)
        x == 1 && return one(T)
        return HUGE_VAL
    end
    x > (a+1)/(a+b+2) && return 1 - ibetar_st(1 - x, b, a; HUGE_VAL, FDeps, rtol, MAXITERS)
    x ≤ FDeps && return zero(T)
    p1 = zero(T)
    q1 = one(T)
    p2 = exp(a*log(x) + b*log(1-x) + loggamma(a+b) - loggamma(a) - loggamma(b)) / a
    q2 = one(T)
    previous = HUGE_VAL
    k = 0
    while k < MAXITERS && !isapprox(p2, previous; rtol)
        previous = p2
        d = -(a+k)*(a+b+k)*x / ((a+2k)*(a+2k+1))
        p1 = p1*d + p2
        q1 = q1*d + q2
        k += 1
        d = k*(b-k)*x /((a+2k-1)*(a+2k))
        p2 = p2*d + p1
        q2 = q2*d + q1
        if q2 ≤ FDeps
            p2 = HUGE_VAL
        else
            p1 /= q2
            q1 /= q2
            p2 /= q2
            q2 = one(T)
        end
    end
    !isapprox(p2, previous; rtol) && error("ibetar did not converge.")
    p2
end

# %%
a, b = 100.0, 0.01
x = 0.1
@show incbeta_I(x, a, b)
@show ibetar_st(x, a, b);

# %%
@btime incbeta_I($x, $a, $b)
@btime ibetar_st($x, $a, $b);

# %%
@which SpecialFunctions._beta_inc(a, b, x)

# %%
@code_warntype SpecialFunctions._beta_inc(a, b, x, 1-x)

# %%
epps = max(eps(), 1.0e-15)

# %%
@btime SpecialFunctions.beta_inc_power_series($a, $b, $x, $epps)

# %%
y = 1 - x
lambda = a > b ? (a + b)*y - b : a - (a + b)*x

# %%
SpecialFunctions.beta_inc_asymptotic_symmetric(b, a, lambda, 100.0*eps())

# %%
