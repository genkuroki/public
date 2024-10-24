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

res = O.chisq_test_for_risk_ratio([5 1; 1 4])

display(res)
