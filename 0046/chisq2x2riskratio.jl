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
#     display_name: Julia 1.10.1
#     language: julia
#     name: julia-1.10
# ---

# %%
module O

function chisq_test_for_risk_ratio(data; null=1, siglev=0.05)
    pvalue = pvalue_rr_pearson_chisq(data'...; ρ=null)
    riskratiohat = _riskratiohat(data'...)
    confint = confint_rr_pearson_chisq(data'...; α=siglev)
    
    print("Pearson's chi-squared test for risk ratio\n")
    print("  data: ", data, "\n")
    print("  null hypothesis: risk ratio = ", null, "\n")
    print("  P-value: ", round(pvalue; sigdigits=4), "\n")
    print("  point estimate of risk ratio: ",
        round(riskratiohat; sigdigits=4), "\n")
    print("  $(100(1 - siglev))% confidence interval of risk ratio: ",
        round.(confint; sigdigits=4), "\n")
end

using Distributions
using Roots
using StatsFuns

safemul(x, y) = x==0 ? zero(x/y) : x*y
safediv(x, y) = x==0 ? zero(x/y) : x/y

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

function pvalue_rr_pearson_chisq(a, b, c, d; ρ=1.0)
    χ² = chisqstat_rr(a, b, c, d; ρ)
    ccdf(Chisq(1), χ²)
end

function confint_rr_pearson_chisq(a, b, c, d; α=0.05)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0.0, Inf]
    f(logρ) = logit(pvalue_rr_pearson_chisq(a, b, c, d; ρ=exp(logρ))) - logit(α)
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
    [exp(L), exp(U)]
end

end

# %%
O.chisq_test_for_risk_ratio([5 1; 1 4])

# %%
module O

@kwdef struct Chisq2x2RiskRatio{D, N, P, R, L, I}
    data::D
    null::N
    pvalue::P
    riskratiohat::R
    siglev::L
    confint::I
end

function chisq_test_for_risk_ratio(data; null=1, siglev=0.05)
    pvalue = pvalue_rr_pearson_chisq(data'...; ρ=null)
    riskratiohat = _riskratiohat(data'...)
    confint = confint_rr_pearson_chisq(data'...; α=siglev)
    Chisq2x2RiskRatio(data, null, pvalue, riskratiohat, siglev, confint)
end

function Base.show(io::IO, ::MIME"text/plain", x::Chisq2x2RiskRatio)
    (; data, null, pvalue, riskratiohat, siglev, confint) = x
    print(io, "Pearson's chi-squared test for risk ratio\n")
    print(io, "  data: ", data, "\n")
    print(io, "  null hypothesis: risk ratio = ", null, "\n")
    print(io, "  P-value: ", round(pvalue; sigdigits=4), "\n")
    print(io, "  point estimate of risk ratio: ",
        round(riskratiohat; sigdigits=4), "\n")
    print(io, "  $(100(1 - siglev))% confidence interval of risk ratio: ",
        round.(confint; sigdigits=4), "\n")
end

Base.:(==)(x::Chisq2x2RiskRatio, y::Chisq2x2RiskRatio) =
    all(fieldnames(Chisq2x2RiskRatio)) do k
        getproperty(x, k) == getproperty(y, k)
    end

function Base.show(io::IO, x::Chisq2x2RiskRatio)
    (; data, null, pvalue, riskratiohat, siglev, confint) = x
    print(io,
        "Chisq2x2RiskRatio(",
        "data=", data, ", ",
        "null=", null, ", ",
        "pvalue=", pvalue, ", ",
        "riskratiohat=", riskratiohat, ", ",
        "siglev=", siglev, ", ",
        "confint=", confint, ")")
end

using Distributions
using Roots
using StatsFuns

safemul(x, y) = x==0 ? zero(x/y) : x*y
safediv(x, y) = x==0 ? zero(x/y) : x/y

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

function pvalue_rr_pearson_chisq(a, b, c, d; ρ=1.0)
    χ² = chisqstat_rr(a, b, c, d; ρ)
    ccdf(Chisq(1), χ²)
end

function confint_rr_pearson_chisq(a, b, c, d; α=0.05)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0.0, Inf]
    f(logρ) = logit(pvalue_rr_pearson_chisq(a, b, c, d; ρ=exp(logρ))) - logit(α)
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
    [exp(L), exp(U)]
end

end

# %%
res = O.chisq_test_for_risk_ratio([5 1; 1 4])

# %%
dump(res)

# %%
show(res)

# %%
display(res)

# %%
string(res)

# %%
res2 = eval(Meta.parse("O." * string(res)))

# %%
res == res2

# %%
foo_jl = read("foo.jl", String)
display(MIME"text/markdown"(), "```julia\n$(foo_jl)```")

# %%
; julia foo.jl

# %%
using StatsPlots

# %%
A = [5 1; 1 4]
@show log₁₀ci = log10.(O.confint_rr_pearson_chisq(A'...; α=0.05))
plot(x -> O.pvalue_rr_pearson_chisq(A'...; ρ=10^x), -1, 2;
    label="", ylim=(-0.03, 1.03), xguide="log₁₀(risk ratio)", yguide="P-value",
    ytick=0:0.05:1) |> display
res = O.chisq_test_for_risk_ratio(A)

# %%
A = [3 0; 0 3]
@show log₁₀ci = log10.(O.confint_rr_pearson_chisq(A'...; α=0.05))
plot(x -> O.pvalue_rr_pearson_chisq(A'...; ρ=10^x), -2, 10;
    label="", ylim=(-0.03, 1.03), xguide="log₁₀(risk ratio)", yguide="P-value",
    ytick=0:0.05:1) |> display
res = O.chisq_test_for_risk_ratio(A)

# %%
A = [3 0; 0.001 3]
@show log₁₀ci = log10.(O.confint_rr_pearson_chisq(A'...; α=0.05))
plot(x -> O.pvalue_rr_pearson_chisq(A'...; ρ=10^x), -2, 10;
    label="", ylim=(-0.03, 1.03), xguide="log₁₀(risk ratio)", yguide="P-value",
    ytick=0:0.05:1) |> display
res = O.chisq_test_for_risk_ratio(A)

# %%
A = [0 1; 1 1]
@show log₁₀ci = log10.(O.confint_rr_pearson_chisq(A'...; α=0.05))
plot(x -> O.pvalue_rr_pearson_chisq(A'...; ρ=10^x), -10, 2;
    label="", ylim=(-0.03, 1.03), xguide="log₁₀(risk ratio)", yguide="P-value",
    ytick=0:0.05:1) |> display
res = O.chisq_test_for_risk_ratio(A)

# %%
A = [0.0001 1; 1 1]
@show log₁₀ci = log10.(O.confint_rr_pearson_chisq(A'...; α=0.05))
plot(x -> O.pvalue_rr_pearson_chisq(A'...; ρ=10^x), -10, 2;
    label="", ylim=(-0.03, 1.03), xguide="log₁₀(risk ratio)", yguide="P-value",
    ytick=0:0.05:1) |> display
res = O.chisq_test_for_risk_ratio(A)

# %%
A = [1 1; 0 1]
@show log₁₀ci = log10.(O.confint_rr_pearson_chisq(A'...; α=0.05))
plot(x -> O.pvalue_rr_pearson_chisq(A'...; ρ=10^x), -2, 10;
    label="", ylim=(-0.03, 1.03), xguide="log₁₀(risk ratio)", yguide="P-value",
    ytick=0:0.05:1) |> display
res = O.chisq_test_for_risk_ratio(A)

# %%
A = [1 1; 0.0001 1]
@show log₁₀ci = log10.(O.confint_rr_pearson_chisq(A'...; α=0.05))
plot(x -> O.pvalue_rr_pearson_chisq(A'...; ρ=10^x), -2, 10;
    label="", ylim=(-0.03, 1.03), xguide="log₁₀(risk ratio)", yguide="P-value",
    ytick=0:0.05:1) |> display
res = O.chisq_test_for_risk_ratio(A)

# %%
