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


* Gilmar Reis, M.D., Ph.D., Eduardo A.S.M. Silva, M.D., Ph.D., Daniela C.M. Silva, M.D., Ph.D., Lehana Thabane, Ph.D., Aline C. Milagres, R.N., Thiago S. Ferreira, M.D., Castilho V.Q. dos Santos, Vitoria H.S. Campos, Ana M.R. Nogueira, M.D., Ana P.F.G. de Almeida, M.D., Eduardo D. Callegari, M.D., Adhemar D.F. Neto, M.D., Ph.D., et al., for the TOGETHER Investigators. Effect of Early Treatment with Ivermectin among Patients with Covid-19. \[[html](https://www.nejm.org/doi/full/10.1056/NEJMoa2115869)\]
<!-- #endregion -->

```julia
using Distributions
using StatsPlots
default(fmt = :png, titlefontsize=12)
using KernelDensity
using Printf
using Roots
```

## Supplementary Appendix Figure S6 の再現

* [Supplementary Appendix](https://www.nejm.org/doi/suppl/10.1056/NEJMoa2115869/suppl_file/nejmoa2115869_appendix.pdf)

のp.15のFigure S6の再現

```julia
# Ivermectin  a b
# Placebo     c d

function bayesian_binomial(A;
        α = 1.0, β = 1.0,
        nsims = 10^6,
        alpha = 0.05,
        title = "",
        ERlim = (0.1, 0.25),
        RRlim = (0.5, 2.5),
        ERtick = 0:0.025:1,
        RRtick = 0:0.25:10,
    )
    a, b, c, d = A'
    beta1 = Beta(α + a, β + b)
    beta2 = Beta(α + c, β + d)
    
    R1 = rand(beta1, nsims)
    R2 = rand(beta2, nsims)
    RR = R1 ./ R2
    L, M, U = quantile.(Ref(RR), (alpha/2, 0.5, 1 - alpha/2))
    Lstr = @sprintf "%.2f" L
    Mstr = @sprintf "%.2f" M
    Ustr = @sprintf "%.2f" U
    BCIstr = "RR [BCI]: $Mstr [$Lstr, $Ustr]"
    ik = InterpKDE(kde(RR))
    f(x) = pdf(ik, x)
    
    P1 = plot(; title)
    plot!(x -> pdf(beta1, x), ERlim...; c=:red, label="Ivermectin")
    plot!(x -> pdf(beta2, x), ERlim...; c=:blue, label="Placebo")
    plot!(; xlabel="Event rate", ylabel="Posterior Density")
    plot!(; xtick = ERtick)
    
    P2 = plot(; title = "$BCIstr")
    plot!(f, L, U; c=:black, lw=2, label="", fillrange=0, fillcolor=:cyan)
    plot!(f, first(RRlim), L; c=:black, lw=2, label="", fillrange=0, fillcolor=:red)
    plot!(f, U, last(RRlim);  c=:black, lw=2, label="", fillrange=0, fillcolor=:red)
    plot!(; xlabel="Relative risk", ylabel="Posterior Density")
    plot!(; xtick = RRtick)
    
    plot(P1, P2; size=(1000, 300))
    plot!(; leftmargin=5Plots.mm, bottommargin=5Plots.mm)
end
```

```julia
bayesian_binomial([100 679-100; 111 679-111];
    title = "Intention-to-treat analysis",
    ERlim = (0.105, 0.215),
    RRlim = (0.50, 1.70)
) |> display

bayesian_binomial([95 674-95; 107 675-107];
    title = "Modified intention-to-treat analysis",
    ERlim = (0.105, 0.21),
    RRlim = (0.50, 1.70)
) |> display

bayesian_binomial([82 624-82; 40 288-40];
    title = "Per-protocol analysis",
    ERlim = (0.08, 0.22),
    RRlim = (0.40, 2.60),
    ERtick = 0:0.03:1
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

function Delta(a, b, c, d, ρ)
    m, n = a + b, c + d
    A, B, C = 1 - ρ, n - a + ρ*(m - c), -a*n + ρ*m*c
    safediv(2C, -B - safesqrt(B^2 - 4A*C))
end
Delta(A, ρ) = Delta(A'..., ρ)

function varDelta(a, b, c, d, ρ)
    Δ = Delta(a, b, c, d, ρ)
    1/(1/(a - Δ) - 1/(a + b - Δ) + 1/(c + Δ) - 1/(c + d + Δ))
end
varDelta(A::AbstractVecOrMat, ρ) = varDelta(A'..., ρ)

chisqstat(A, ρ = 1.0) = safediv(Delta(A, ρ)^2, varDelta(A, ρ))
pvalue_chisq(A, ρ = 1.0) = ccdf(Chisq(1), chisqstat(A, ρ))
function ci_chisq(A, α = 0.05)
    f(t) = pvalue_chisq(A, exp(t)) - α
    logci = find_zeros(f, -1e1, 1e1)
    exp(first(logci)), exp(last(logci))
end

function pvalue_bin(a, b, p)
    bin = Binomial(a + b, p)
    min(1, 2cdf(bin, a), 2ccdf(bin, a-1))
end

function chisq_test(A; RR₀ = 1.0, alpha = 0.05)
    RR = riskratio(A)
    chisq = chisqstat(A, RR₀)
    df = 1
    p_value = ccdf(Chisq(df), chisq)
    conf_int = ci_chisq(A, alpha)
    (; RR, RR₀, p_value, alpha, conf_int, chisq, df)
end

function plot_chisq_test(A; RR₀ = 1.0, alpha = 0.05,
        title = "",
        ERlim = (0.1, 0.25),
        RRlim = (0.5, 2.5),
        ERtick = 0:0.025:1,
        RRtick = 0:0.25:10,
    )
    (; RR, RR₀, p_value, alpha, conf_int, chisq, df) = chisq_test(A; RR₀, alpha)
    a, b, c, d = A'
    
    f(p) = pvalue_bin(a, b, p)
    g(q) = pvalue_bin(c, d, q)
    
    h(RR) = pvalue_chisq(A, RR)
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
plot_chisq_test([100 679-100; 111 679-111];
    title = "Intention-to-treat analysis",
    ERlim = (0.105, 0.215),
    RRlim = (0.50, 1.70)
) |> display

plot_chisq_test([95 674-95; 107 675-107];
    title = "Modified intention-to-treat analysis",
    ERlim = (0.105, 0.21),
    RRlim = (0.50, 1.70)
) |> display

plot_chisq_test([82 624-82; 40 288-40];
    title = "Per-protocol analysis",
    ERlim = (0.08, 0.22),
    RRlim = (0.40, 2.60),
    ERtick = 0:0.03:1
) |> display
```

```julia

```
