---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.3
  kernelspec:
    display_name: Julia 1.8.0-beta1
    language: julia
    name: julia-1.8
---

<!-- #region -->
# イベルメクチン論文の図の再現

* 黒木玄
* 2022-04-02


* Gilmar Reis, M.D., Ph.D., Eduardo A.S.M. Silva, M.D., Ph.D., Daniela C.M. Silva, M.D., Ph.D., Lehana Thabane, Ph.D., Aline C. Milagres, R.N., Thiago S. Ferreira, M.D., Castilho V.Q. dos Santos, Vitoria H.S. Campos, Ana M.R. Nogueira, M.D., Ana P.F.G. de Almeida, M.D., Eduardo D. Callegari, M.D., Adhemar D.F. Neto, M.D., Ph.D., et al., for the TOGETHER Investigators. Effect of Early Treatment with Ivermectin among Patients with Covid-19. New England Journal of Medicine, March 30, 2022
DOI: 10.1056/NEJMoa2115869 \[[html](https://www.nejm.org/doi/full/10.1056/NEJMoa2115869)\] \[[pdf](https://www.nejm.org/doi/pdf/10.1056/NEJMoa2115869?articleTools=true)\]
<!-- #endregion -->

```julia
using Distributions
using StatsPlots
default(fmt = :png, titlefontsize=12)
using KernelDensity
using Printf
using Roots
using QuadGK
```

## Supplementary Appendix Figure S6 の再現

* Supplementary Appendix \[[pdf](https://www.nejm.org/doi/suppl/10.1056/NEJMoa2115869/suppl_file/nejmoa2115869_appendix.pdf)\]

のp.15のFigure S6の再現. データは[論文](https://www.nejm.org/doi/full/10.1056/NEJMoa2115869)のTable 2より.

```julia
# Ivermectin  a b
# Placebo     c d

function pdfRR(beta1, beta2, ρ)
    f(q) = exp(logpdf(beta1, ρ*q) + log(q) + logpdf(beta2, q)) 
    if 0 < ρ ≤ 1
        quadgk(f, 0, 1)[1]
    elseif ρ > 1
        quadgk(f, 0, 1/ρ)[1]
    else
        zero(ρ)
    end
end

function cdfRR(beta1, beta2, ρ)
    if 0 < ρ ≤ 1
        f(p) = exp(logpdf(beta1, p) + logccdf(beta2, p/ρ))
        quadgk(f, 0, ρ)[1]
    elseif ρ > 1
        g(q) = exp(logpdf(beta2, q) + logccdf(beta1, ρ*q))
        1 - quadgk(g, 0, 1/ρ)[1]
    else
        zero(ρ)
    end
end

function bayesian_binomial(A;
        α = 1.0, β = 1.0,
        alpha = 0.05,
        nsims = 10^6,
        title = "",
        ERlim = (0.1, 0.25),
        RRlim = (0.5, 2.5),
        ERtick = 0:0.025:1,
        RRtick = 0:0.25:10,
    )
    a, b, c, d = A'
    beta1 = Beta(α + a, β + b)
    beta2 = Beta(α + c, β + d)
    
    if nsims == 0
        lw = 1
        L = find_zero(ρ -> cdfRR(beta1, beta2, ρ) - alpha/2, 1.0)
        M = find_zero(ρ -> cdfRR(beta1, beta2, ρ) - 0.5, 1.0)
        U = find_zero(ρ -> cdfRR(beta1, beta2, ρ) - (1 - alpha/2), 1.0)
        f = ρ -> pdfRR(beta1, beta2, ρ)
    else
        lw = 2
        R1, R2 = rand(beta1, nsims), rand(beta2, nsims)
        RR = R1 ./ R2
        L, M, U = quantile.(Ref(RR), (alpha/2, 0.5, 1 - alpha/2))
        ik = InterpKDE(kde(RR))
        f = ρ -> pdf(ik, ρ)
    end
    Lstr = @sprintf "%.2f" L
    Mstr = @sprintf "%.2f" M
    Ustr = @sprintf "%.2f" U
    BCIstr = "RR [BCI]: $Mstr [$Lstr, $Ustr]"
     
    P1 = plot(; title)
    plot!(x -> pdf(beta1, x), ERlim...; c=:red, label="Ivermectin")
    plot!(x -> pdf(beta2, x), ERlim...; c=:blue, label="Placebo")
    plot!(; xlabel="Event rate", ylabel="Posterior Density")
    plot!(; xtick = ERtick)
    
    P2 = plot(; title = "$BCIstr")
    plot!(f, L, U; c=:black, lw, label="", fillrange=0, fillcolor=:cyan)
    plot!(f, first(RRlim), L; c=:black, lw, label="", fillrange=0, fillcolor=:red)
    plot!(f, U, last(RRlim);  c=:black, lw, label="", fillrange=0, fillcolor=:red)
    plot!(; xlabel="Relative risk", ylabel="Posterior Density")
    plot!(; xtick = RRtick)
    
    plot(P1, P2; size=(1000, 300))
    plot!(; leftmargin=5Plots.mm, bottommargin=5Plots.mm)
end
```

```julia
nsims = 0

bayesian_binomial([100 679-100; 111 679-111]; nsims,
    title = "Intention-to-treat analysis",
    ERlim = (0.105, 0.215), RRlim = (0.50, 1.70),
) |> display

bayesian_binomial([95 674-95; 107 675-107]; nsims,
    title = "Modified intention-to-treat analysis",
    ERlim = (0.105, 0.21), RRlim = (0.50, 1.70),
) |> display

bayesian_binomial([82 624-82; 40 288-40]; nsims,
    title = "Per-protocol analysis",
    ERlim = (0.08, 0.22), RRlim = (0.40, 2.60), ERtick = 0:0.03:1,
) |> display
```

```julia
nsims = 10^6

bayesian_binomial([100 679-100; 111 679-111]; nsims,
    title = "Intention-to-treat analysis",
    ERlim = (0.105, 0.215), RRlim = (0.50, 1.70),
) |> display

bayesian_binomial([95 674-95; 107 675-107]; nsims,
    title = "Modified intention-to-treat analysis",
    ERlim = (0.105, 0.21), RRlim = (0.50, 1.70),
) |> display

bayesian_binomial([82 624-82; 40 288-40]; nsims,
    title = "Per-protocol analysis",
    ERlim = (0.08, 0.22), RRlim = (0.40, 2.60), ERtick = 0:0.03:1,
) |> display
```

## P値版のグラフ

```julia
# Ivermectin  a b
# Placebo     c d

safediv(x, y) = x == 0 ? x/one(y) : x/y
safesqrt(x) = √max(0, x)

riskratio(a, b, c, d) = safediv(a*(c+d), (a+b)*c)
riskratio(A) = riskratio(A'...)

"""
Δ = Delta(a, b, c, d, ρ) is characterized by ((a - Δ)(c + Δ + d))/((a - Δ + b)(c + Δ)) = ρ.
"""
function Delta(a, b, c, d, ρ)
    m, n = a + b, c + d
    A, B, C = 1 - ρ, n - a + ρ*(m - c), -a*n + ρ*m*c
    safediv(2C, -B - safesqrt(B^2 - 4A*C))
end
Delta(A, ρ) = Delta(A'..., ρ)

"""estimator of var(Δ)"""
function varDelta(a, b, c, d, ρ)
    Δ = Delta(a, b, c, d, ρ)
    1/(1/(a - Δ) - 1/(a + b - Δ) + 1/(c + Δ) - 1/(c + d + Δ))
end
varDelta(A::AbstractVecOrMat, ρ) = varDelta(A'..., ρ)

chisq_RR(A, ρ = 1.0) = safediv(Delta(A, ρ)^2, varDelta(A, ρ))

"""chi-squared P-value function of relatice risk"""
pvalue_chisq_RR(A, ρ = 1.0) = ccdf(Chisq(1), chisq_RR(A, ρ))

"""chi-squared confidence interval of relatice risk"""
function ci_chisq_RR(A, α = 0.05)
    f(t) = pvalue_chisq_RR(A, exp(t)) - α
    logci = find_zeros(f, -1e1, 1e1)
    exp(first(logci)), exp(last(logci))
end

"""Clopper-Pearson P-value function of binomial distribution model"""
function pvalue_clopper_pearson(a, b, p)
    bin = Binomial(a + b, p)
    min(1, 2cdf(bin, a), 2ccdf(bin, a-1))
end

"""chi-squared P-value function of binomial distribution model"""
function pvalue_chisq_bin(a, b, p)
    m = a + b
    chisq = safediv((a - m*p)^2, m*p*(1 - p))
    ccdf(Chisq(1), chisq)
end

function chisq_test_RR(A; RR₀ = 1.0, alpha = 0.05)
    RR = riskratio(A)
    chisq = chisq_RR(A, RR₀)
    df = 1
    p_value = ccdf(Chisq(df), chisq)
    conf_int = ci_chisq_RR(A, alpha)
    (; RR, RR₀, p_value, alpha, conf_int, chisq, df)
end

function plot_chisq_test(A; RR₀ = 1.0, alpha = 0.05,
        title = "",
        ERlim = (0.1, 0.25),
        RRlim = (0.5, 2.5),
        ERtick = 0:0.025:1,
        RRtick = 0:0.25:10,
        pvalue_bin = pvalue_chisq_bin,
    )
    (; RR, RR₀, p_value, alpha, conf_int, chisq, df) = chisq_test_RR(A; RR₀, alpha)
    a, b, c, d = A'
    
    f(p) = pvalue_bin(a, b, p)
    g(q) = pvalue_bin(c, d, q)
    
    h(RR) = pvalue_chisq_RR(A, RR)
    L, U = conf_int
    Lstr = @sprintf "%.2f" L
    RRstr = @sprintf "%.2f" RR
    Ustr = @sprintf "%.2f" U
    CIstr = "RR [CI]: $RRstr [$Lstr, $Ustr]"
    
    P1 = plot(; title)
    plot!(f, ERlim...; c=:red, label="Ivermectin")
    plot!(g, ERlim...; c=:blue, label="Placebo")
    plot!(; xlabel="Event rate", ylabel="P-value")
    plot!(; xtick = ERtick, ytick=0:0.1:1)
    
    P2 = plot(; title = "$CIstr")
    plot!(h, RRlim...; c=:black, label="")
    plot!([L, U], [alpha, alpha]; c=:red, lw=4, label="")
    plot!(; xlabel="Relative risk", ylabel="P-value")
    plot!(; xtick = RRtick, ytick=0:0.1:1)
    
    plot(P1, P2; size=(1000, 300))
    plot!(; leftmargin=5Plots.mm, bottommargin=5Plots.mm)
end

```

```julia
pvalue_bin = pvalue_chisq_bin

plot_chisq_test([100 679-100; 111 679-111]; pvalue_bin,
    title = "Intention-to-treat analysis",
    ERlim = (0.105, 0.215), RRlim = (0.50, 1.70),
) |> display

plot_chisq_test([95 674-95; 107 675-107]; pvalue_bin,
    title = "Modified intention-to-treat analysis",
    ERlim = (0.105, 0.21), RRlim = (0.50, 1.70),
) |> display

plot_chisq_test([82 624-82; 40 288-40]; pvalue_bin,
    title = "Per-protocol analysis",
    ERlim = (0.08, 0.22), RRlim = (0.40, 2.60), ERtick = 0:0.03:1,
) |> display
```

```julia
pvalue_bin = pvalue_clopper_pearson

plot_chisq_test([100 679-100; 111 679-111]; pvalue_bin,
    title = "Intention-to-treat analysis",
    ERlim = (0.105, 0.215), RRlim = (0.50, 1.70),
) |> display

plot_chisq_test([95 674-95; 107 675-107]; pvalue_bin,
    title = "Modified intention-to-treat analysis",
    ERlim = (0.105, 0.21), RRlim = (0.50, 1.70),
) |> display

plot_chisq_test([82 624-82; 40 288-40]; pvalue_bin,
    title = "Per-protocol analysis",
    ERlim = (0.08, 0.22), RRlim = (0.40, 2.60), ERtick = 0:0.03:1,
) |> display
```

```julia

```
