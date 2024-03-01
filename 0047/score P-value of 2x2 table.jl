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
using Distributions
using Roots
using StatsPlots
default(fmt=:png)

myecdf(A, x) = count(≤(x), A)/length(A)
safemul(x, y) = x == 0 ? zero(x*y) : y == 0 ? zero(x*y) : x*y
safediv(x, y) = x == 0 ? zero(x/y) : isinf(y) ? zero(x/y) : x/y

# %%
### score P-value for rate difference

riskdiffhat_score(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function loglik_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safemul(a, log(p)) + safemul(b, log(1-p)) + safemul(c, log(q)) + safemul(d, log(1-q))
end

function scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(a, p) - safediv(b, 1-p) + safediv(c, q) - safediv(d, 1-q)
end

function d_scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    -safediv(a, p^2) - safediv(b, (1-p)^2) - safediv(c, q^2) - safediv(d, (1-q)^2)
end

function scorestat_Δ_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(a, p) - safediv(b, 1-p)
end

function estimate_q_given_Δ_rd(a, b, c, d, Δ=0.0; alg=Bisection())
    qmin, qmax = max(0.0, -Δ), min(1.0, 1.0-Δ)
    a+c==0 && return qmin
    b+d==0 && return qmax
    f(q) = scorestat_q_rd(a, b, c, d, q, Δ)
    S_qmin = f(qmin + eps())
    S_qmax = f(qmax - eps())
    S_qmin ≥ 0 && S_qmax ≥ 0 && return S_qmin < S_qmax ? qmin : qmax
    S_qmin ≤ 0 && S_qmax ≤ 0 && return S_qmin < S_qmax ? qmax : qmin
    find_zero(f, (qmin + eps(), qmax - eps()), alg)
end

function varinv_scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(p*(1-p), a+b) + safediv(q*(1-q), c+d)
end

function chisqstat_rd_score(a, b, c, d; Δ=0.0, alg=Bisection())
    abs(Δ) == 1 && return Inf
    q̃ = estimate_q_given_Δ_rd(a, b, c, d, Δ; alg)
    S = scorestat_Δ_rd(a, b, c, d, q̃, Δ)
    Vinv = varinv_scorestat_q_rd(a, b, c, d, q̃, Δ)
    safemul(S^2, Vinv)
end

function pvalue_rd_score(a, b, c, d; Δ=0.0, alg=Bisection())
    χ² = chisqstat_rd_score(a, b, c, d; Δ, alg)
    ccdf(Chisq(1), χ²)
end

function confint_rd_score(a, b, c, d; α=0.05, alg=Bisection())
    χ²_α = cquantile(Chisq(1), α)
    g(Δ) = chisqstat_rd_score(a, b, c, d; Δ, alg) - χ²_α
    Δ0 = riskdiffhat_score(a, b, c, d)
    L = if g(-1 + eps()) > 0
        find_zero(g, (-1 + eps(), Δ0), alg)
    else
        -1.0
    end
    U = if g(1 - eps()) > 0
        find_zero(g, (Δ0, 1 - eps()), alg)
    else
        1.0
    end
    [L, U]
end

# %%
a, b, c, d = 8, 0, 8, 0
@show Δ0 = riskdiffhat_score(a, b, c, d)
@show chisqstat_rd_score(a, b, c, d; Δ=Δ0)
chisqstat_rd_score(a, b, c, d; Δ=-0.9999)

# %%
a, b, c, d = 1, 2, 3, 4
pvalue_rd_score(1, 2, 3, 4), ccdf(Chisq(1), (a+b+c+d)*(a*d-b*c)^2/((a+b)*(c+d)*(a+c)*(b+d)))

# %%
a, b, c, d = 3, 0, 4, 0
pvalue_rd_score(a, b, c, d), ccdf(Chisq(1), safediv((a+b+c+d)*(a*d-b*c)^2, (a+b)*(c+d)*(a+c)*(b+d)))

# %%
riskdiffhat(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function stderr_riskdiffhat(a, b, c, d)
    m, n = a+b, c+d
    p̂, q̂ = safediv(a, m), safediv(c, n)
    √(safediv(p̂*(1-p̂), m) + safediv(q̂*(1-q̂), n))
end

function pvalue_rd_wald(a, b, c, d; Δ=0)
    RDhat = riskdiffhat(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat(a, b, c, d)
    2ccdf(Normal(0, 1), safediv(abs(RDhat - Δ), SEhat_riskdiffhat))
end

function confint_rd_wald(a, b, c, d; α=0.05)
    z = quantile(Normal(), 1-α/2)
    RDhat = riskdiffhat(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat(a, b, c, d)
    [RDhat - z*SEhat_riskdiffhat, RDhat + z*SEhat_riskdiffhat]
end

riskdiffhat_zou_donner(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function stderr_riskdiffhat_zou_donner(a, b, c, d; u=1)
    m, n = a+b, c+d
    p̂, q̂ = safediv(a, m), safediv(c, n)
    √(safediv(p̂*(1-p̂), m-u) + safediv(q̂*(1-q̂), n-u))
end

function pvalue_rd_zou_donner(a, b, c, d; Δ=0, u=1)
    ((a==0 && d==0) || (b==0 && c==0)) && return 1.0
    RDhat = riskdiffhat_zou_donner(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat_zou_donner(a, b, c, d; u)
    Z = safediv((1 - RDhat^2)*abs(atanh(RDhat) - atanh(Δ)), SEhat_riskdiffhat)
    2ccdf(Normal(), abs(Z))
end

function confint_rd_zou_donner(a, b, c, d; α=0.05, u=1)
    z = quantile(Normal(), 1-α/2)
    RDhat = riskdiffhat_zou_donner(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat_zou_donner(a, b, c, d; u)
    m = atanh(RDhat)
    d = safediv(z*SEhat_riskdiffhat, 1 - RDhat^2)
    [tanh(m-d), tanh(m+d)]
end

a, b = 58, 22
c, d = 62, 38

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show confint_rd_wald(a, b, c, d; α=0.05)
@show confint_rd_zou_donner(a, b, c, d; α=0.05)
@show confint_rd_score(a, b, c, d; α=0.05)

plot(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -0.2, 0.4; label="")
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ), -0.2, 0.4; label="", ls=:dash)
plot!(Δ -> pvalue_rd_score(a, b, c, d; Δ), -0.2, 0.4; label="", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-2:0.05:2, ytick=0:0.05:1)

# %%
a, b = 8, 2
c, d = 3, 7

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=-1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=-1.0)
@show pvalue_rd_score(a, b, c, d; Δ=-1.0)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show pvalue_rd_wald(a, b, c, d; Δ=1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0)
@show pvalue_rd_score(a, b, c, d; Δ=1.0)
@show confint_rd_wald(a, b, c, d; α=0.05)
@show confint_rd_zou_donner(a, b, c, d; α=0.05)
@show confint_rd_score(a, b, c, d; α=0.05)

plot(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -1, 1.0; label="Wald")
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ); label="Zou-Donner", ls=:dash)
plot!(Δ -> pvalue_rd_score(a, b, c, d; Δ); label="score", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-2:0.1:2, ytick=0:0.05:1)

# %%
a, b = 190, 10
c, d = 180, 20

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=-1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=-1.0)
@show pvalue_rd_score(a, b, c, d; Δ=-1.0)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show pvalue_rd_wald(a, b, c, d; Δ=1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0)
@show pvalue_rd_score(a, b, c, d; Δ=1.0)
@show confint_rd_wald(a, b, c, d; α=0.05)
@show confint_rd_zou_donner(a, b, c, d; α=0.05)
@show confint_rd_score(a, b, c, d; α=0.05)

plot(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -0.5, 1.0; label="Wald")
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ); label="Zou-Donner", ls=:dash)
plot!(Δ -> pvalue_rd_score(a, b, c, d; Δ); label="score", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-2:0.1:2, ytick=0:0.05:1)

# %%
a, b = 8, 0
c, d = 9, 0

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=-1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=-1.0)
@show pvalue_rd_score(a, b, c, d; Δ=-1.0)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show pvalue_rd_wald(a, b, c, d; Δ=1.0-eps())
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0-eps())
@show pvalue_rd_score(a, b, c, d; Δ=1.0-eps())
@show pvalue_rd_wald(a, b, c, d; Δ=1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0)
@show pvalue_rd_score(a, b, c, d; Δ=1.0)
@show confint_rd_wald(a, b, c, d; α=0.05)
@show confint_rd_zou_donner(a, b, c, d; α=0.05)
@show confint_rd_score(a, b, c, d; α=0.05)

plot(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -1.0, 1.0; label="Wald")
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ); label="Zou-Donner", ls=:dash)
plot!(Δ -> pvalue_rd_score(a, b, c, d; Δ); label="score", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-2:0.1:2, ytick=0:0.05:1)

# %%
a, b = 190, 10
c, d = 180, 20
Δ = 0.5

@show estimate_q_given_Δ_rd(a, b, c, d, Δ)
@show scorestat_q_rd(a, b, c, d, 0.0, Δ)
@show scorestat_q_rd(a, b, c, d, 0.675, Δ)
@show scorestat_q_rd(a, b, c, d, 0.7, Δ)

find_zero(q -> (q; scorestat_q_rd(a, b, c, d, q, Δ)), (0.0, 1-Δ))

# %%
a, b = 190, 10
c, d = 180, 20
Δ = -1.0

@show scorestat_q_rd(a, b, c, d, 0.0, 0.3)
@show scorestat_q_rd(a, b, c, d, 0.0, 0.3)

qmin, qmax = max(0.0, -Δ), min(1.0, 1.0-Δ)
a+c==0 && @show qmin
b+d==0 && @show qmax
@show S_qmin = scorestat_q_rd(a, b, c, d, qmin+eps(), Δ)
@show S_qmax = scorestat_q_rd(a, b, c, d, qmax-eps(), Δ)
S_qmin ≥ 0 && S_qmax ≥ 0 && @show S_qmin < S_qmax ? qmin : qmax
S_qmin ≤ 0 && S_qmax ≤ 0 && @show S_qmin < S_qmax ? qmax : qmin
f(q) = scorestat_q_rd(a, b, c, d, q, Δ)
@show q0 = (qmin + qmax)/2
@show find_zero(f, q0, Order0())
@show find_zero(f, (qmin+eps(), qmax-eps()))


plot(q -> scorestat_q_rd(a, b, c, d, q, 0.3), 0.01, 0.69)

# %%
a, b = 10, 0
c, d = 20, 0

@show scorestat_q_rd(a, b, c, d, 0.0, 0.3)
@show scorestat_q_rd(a, b, c, d, 0.7, 0.3)

plot(q -> scorestat_q_rd(a, b, c, d, q, 0.3), 0.1, 0.69)

# %%
a, b = 0, 2
c, d = 0, 7

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)

plot(Δ -> pvalue_rd_score(a, b, c, d; Δ), -1, 1; label="score", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-1:0.2:1, ytick=0:0.05:1)

# %%
a, b = 4, 0
c, d = 10, 1

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show pvalue_rd_wald(a, b, c, d; Δ=-1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=-1.0)
@show pvalue_rd_score(a, b, c, d; Δ=-1.0+√eps())
@show pvalue_rd_wald(a, b, c, d; Δ=1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0)
@show pvalue_rd_score(a, b, c, d; Δ=1.0)

plot(Δ -> pvalue_rd_score(a, b, c, d; Δ), -1, 1; label="score", ls=:dashdot)
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-1:0.2:1, ytick=0:0.05:1)

# %%
a, b = 8, 0
c, d = 8, 1

@show scorestat_q_rd(a, b, c, d, 0.0, 0.2)
@show scorestat_q_rd(a, b, c, d, 0.8, 0.2)

plot(q -> scorestat_q_rd(a, b, c, d, q, 0.2), 0.2, 0.79)

# %%

# %%

# %%

# %%

# %%
m = n = 10
L = 10^4
a = rand(Binomial(m, 0.8), L)
b = m .- a
c = rand(Binomial(n, 0.3), L)
d = n .- c
pval_wald = pvalue_rd_wald.(a, b, c, d; Δ=0.5)
pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ=0.5)
pval_score = pvalue_rd_score.(a, b, c, d; Δ=0.5)

plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
plot!(α -> myecdf(pval_score, α); label="socre")
plot!(identity; label="", c=:gray, ls=:dot)
plot!(size=(400, 400))

# %%
m, n = 20, 20
L = 10^4
a = rand(Binomial(m, 0.8), L)
b = m .- a
c = rand(Binomial(n, 0.3), L)
d = n .- c
pval_wald = pvalue_rd_wald.(a, b, c, d; Δ=0.5)
pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ=0.5)
pval_score = pvalue_rd_score.(a, b, c, d; Δ=0.5)

plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
plot!(α -> myecdf(pval_score, α); label="socre")
plot!(identity; label="", c=:gray, ls=:dot)
plot!(size=(400, 400))

# %%
m, n = 40, 40
L = 10^4
a = rand(Binomial(m, 0.8), L)
b = m .- a
c = rand(Binomial(n, 0.3), L)
d = n .- c
pval_wald = pvalue_rd_wald.(a, b, c, d; Δ=0.5)
pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ=0.5)
pval_score = pvalue_rd_score.(a, b, c, d; Δ=0.5)

plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
plot!(α -> myecdf(pval_score, α); label="score")
plot!(identity; label="", c=:gray, ls=:dot)
plot!(size=(400, 400))

# %%
m, n = 20, 100
L = 10^4
a = rand(Binomial(m, 0.8), L)
b = m .- a
c = rand(Binomial(n, 0.3), L)
d = n .- c
pval_wald = pvalue_rd_wald.(a, b, c, d; Δ=0.5)
pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ=0.5)
pval_score = pvalue_rd_score.(a, b, c, d; Δ=0.5)

plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
plot!(α -> myecdf(pval_score, α); label="score")
plot!(identity; label="", c=:gray, ls=:dot)
plot!(size=(400, 400))

# %%
m, n = 100, 20
L = 10^4
a = rand(Binomial(m, 0.8), L)
b = m .- a
c = rand(Binomial(n, 0.3), L)
d = n .- c
pval_wald = pvalue_rd_wald.(a, b, c, d; Δ=0.5)
pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ=0.5)
pval_score = pvalue_rd_score.(a, b, c, d; Δ=0.5)

plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
plot!(α -> myecdf(pval_score, α); label="score")
plot!(identity; label="", c=:gray, ls=:dot)
plot!(size=(400, 400))

# %%
m, n = 20, 100
p, q = 0.1, 0.2
Δ = p - q
L = 10^4
a = rand(Binomial(m, p), L)
b = m .- a
c = rand(Binomial(n, q), L)
d = n .- c
pval_wald = pvalue_rd_wald.(a, b, c, d; Δ)
pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ)
pval_score = pvalue_rd_score.(a, b, c, d; Δ)

plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
plot!(α -> myecdf(pval_score, α); label="score")
plot!(identity; label="", c=:gray, ls=:dot)
plot!(size=(400, 400))

# %%
m, n = 100, 20
p, q = 0.1, 0.2
Δ = p - q
L = 10^4
a = rand(Binomial(m, p), L)
b = m .- a
c = rand(Binomial(n, q), L)
d = n .- c
pval_wald = pvalue_rd_wald.(a, b, c, d; Δ)
pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ)
pval_score = pvalue_rd_score.(a, b, c, d; Δ)

plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
plot!(α -> myecdf(pval_score, α); label="score")
plot!(identity; label="", c=:gray, ls=:dot)
plot!(size=(400, 400))

# %%
m, n = 60, 60
p, q = 0.1, 0.2
Δ = p - q
L = 10^4
a = rand(Binomial(m, p), L)
b = m .- a
c = rand(Binomial(n, q), L)
d = n .- c
pval_wald = pvalue_rd_wald.(a, b, c, d; Δ)
pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ)
pval_score = pvalue_rd_score.(a, b, c, d; Δ)

plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
plot!(α -> myecdf(pval_score, α); label="score")
plot!(identity; label="", c=:gray, ls=:dot)
plot!(size=(400, 400))

# %%

# %%

# %%

# %%
function polynom_scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    a*(1-p)*q*(1-q) - b*p*q*(1-q) + c*p*(1-p)*(1-q) - d*p*(1-p)*q
end

# %%
a, b, c, d = 8, 1, 8, 1
@show Δ = 0.2
@show riskdiffhat_score(a, b, c, d)
P = plot(q -> polynom_scorestat_q_rd(a, b, c, d, q, Δ), 0, 1-Δ)
Q = plot(q -> loglik_rd(a, b, c, d, q, Δ), 0, 1-Δ)
plot(P, Q; layout=(1, 2), size=(800, 250))

# %%
plot(q -> loglik_rd(a, b, c, d, q, Δ), 1-Δ-0.1, 1-Δ)

# %%
