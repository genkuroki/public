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

# %%
VERSION

# %% [markdown]
# ## Redefinition of constants (structs)
#
# https://julialang.org/blog/2025/10/julia-1.12-highlights/#redefinition_of_constants_structs
#
# Bindings now participate in the "world age" mechanism previously used for methods. This has the effect that constants and structs can be properly redefined. As an example:

# %%
struct Foo{T} a::T end
f(foo::Foo, x) = ((; a) = foo; a*x)
@show foo = Foo(3)
f(foo, 4)

# %%
struct Foo{T} a::T; b::T end
@show foo = Foo(3, 4)
f(foo, 4)

# %%
f(foo::Foo, x) = ((; a, b) = foo; a*x + b)
f(foo, 4)

# %% [markdown]
# ## Score test of risk ratio for 2x2 tables

# %%
using Distributions
using StatsFuns
using Roots

safemul(x, y) = x == 0 ? zero(x*y) : x*y
safediv(x, y) = x == 0 ? zero(x/y) : x/y

_riskratiohat(a, b, c, d) = safediv(a*(c+d), (a+b)*c)

function Delta(a, b, c, d; ρ=1.0)
    m, n = a+b, c+d
    A, B, C = ρ-1, n-a+ρ*(m-c), a*n-ρ*m*c
    Δ = isinf(ρ) ? oftype(ρ, -c) : ρ==0 ? oftype(ρ, a) : safediv(2C, B + √(B^2 - 4A*C))
end

function _chisqstat_rr(a, b, c, d, Δ)
    m, n = a+b, c+d
    safemul(Δ^2, safediv(b, m*(a-Δ)) + safediv(d, n*(c+Δ)))
end

function chisqstat_rr(a, b, c, d; ρ=1.0)
    Δ = Delta(a, b, c, d; ρ)
    _chisqstat_rr(a, b, c, d, Δ)
end

function pvalue_rr_score(a, b, c, d; ρ=1.0)
    χ² = chisqstat_rr(a, b, c, d; ρ)
    ccdf(Chisq(1), χ²)
end

function confint_rr_score(a, b, c, d; α=0.05)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0.0, Inf]
    f(logρ) = logit(pvalue_rr_score(a, b, c, d; ρ=exp(logρ))) - logit(α)
    L = if f(-Inf) > 0
        -Inf
    else
        logRRhat = log(_riskratiohat(a, b, c, d))
        x0 = logRRhat == -Inf ? -10.0 : logRRhat == Inf ? 10.0 : logRRhat - 1
        find_zero(f, x0)
    end
    U = if f(Inf) > 0
        Inf
    else
        logRRhat = log(_riskratiohat(a, b, c, d))
        x0 = logRRhat == -Inf ? -10.0 : logRRhat == Inf ? 10.0 : logRRhat + 1
        find_zero(f, x0)
    end
    exp(L), exp(U)
end

a, b, c, d = 35, 15, 25, 25
ρ = 1
α = 0.05
@show [a b; c d]
@show ρ
@show α
@show pvalue_rr_score(a, b, c, d; ρ=1)
@show _riskratiohat(a, b, c, d)
@show confint_rr_score(a, b, c, d);

# %%
struct ScoreTest2x2RR{I, R}
    a::I
    b::I
    c::I
    d::I
    ρ::R
    α::R
    pvalue::R
    RRhat::R
    CI::Tuple{R, R}
end

ScoreTest2x2RR(a, b, c, d; ρ=1.0, α=0.05) =
    ScoreTest2x2RR(a, b, c, d, float(ρ), float(α),
        pvalue_rr_score(a, b, c, d; ρ=float(ρ)),
        _riskratiohat(a, b, c, d),
        confint_rr_score(a, b, c, d; α=float(α))
    )

a, b, c, d = 35, 15, 30, 30
ρ = 1
α = 0.05
result = ScoreTest2x2RR(a, b, c, d; ρ, α)

# %%
@show result.a
@show result.b
@show result.c
@show result.d
@show result.ρ
@show result.α
@show result.pvalue
@show result.RRhat
@show result.CI; 

# %%
function Base.show(io::IO, ::MIME"text/plain", x::ScoreTest2x2RR)
    r(x) = round(x; sigdigits=3)
    (; a, b, c, d, ρ, α, pvalue, RRhat, CI) = x
    print(io, "Pearson's chi-squared test for risk ratio\n")
    print(io, "  data: ", [a b; c d], "\n")
    print(io, "  testing risk ratio: ", ρ, "\n")
    print(io, "  P-value: ", r(100pvalue), " %\n")
    print(io, "  point estimate of risk ratio: ", r(RRhat), "\n")
    print(io, "  $(100(1 - α)) % confidence interval of risk ratio: [", r(CI[1]), ", ", r(CI[2]), "]\n")
end

a, b, c, d = 35, 15, 30, 30
ρ = 1
α = 0.05
ScoreTest2x2RR(a, b, c, d; ρ, α)

# %%
function Base.show(io::IO, x::ScoreTest2x2RR)
    r(x) = round(x; sigdigits=3)
    (; a, b, c, d, ρ, α, pvalue, RRhat, CI) = x
    print(io,
        "ScoreTest2x2RR(",
        "a=", a, ", ",
        "b=", b, ", ",
        "c=", c, ", ",
        "d=", d, ", ",
        "ρ=", r(ρ), ", ",
        "α=", r(α), ", ",
        "pvalue=", r(pvalue), ", ",
        "RRhat=", r(RRhat), ", ",
        "CI=", r.(CI), ")")
end

a, b, c, d = 35, 15, 30, 30
ρ = 1
α = 0.05
show(ScoreTest2x2RR(a, b, c, d; ρ, α));

# %% [markdown]
# ## 関数はメソッドの集まり

# %%
double(x) = x + x
double(123)

# %%
double("hoge") # error

# %%
double(x::AbstractString) = x^2
double("hoge")

# %%
double('げ') # error

# %%
double(x::AbstractChar) = double(string(x))
double('げ')

# %%
abstract type AbstractFoo end
struct Foo{T} <: AbstractFoo a::T end
double(Foo(123)) # error

# %%
double(x::AbstractFoo) = Foo(double(x.a))
double(Foo(123))

# %%
Foo("hoge")

# %%
string(Foo("hoge"))

# %%
function Base.show(io::IO, x::AbstractFoo)
    print(io, '⟨')
    show(io, x.a)
    print(io, '⟩')
end
Foo("hoge")

# %%
double(Foo("hoge"))

# %%
string(Foo("hoge"))

# %%
