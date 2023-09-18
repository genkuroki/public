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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using Optim
using Roots
using StatsBase: ecdf
using StatsPlots
default(fmt=:png)

function pvalue_normal_approx(nulldist, x)
    m, s = mean(nulldist), std(nulldist)
    z = (x - m)/s
    2ccdf(Normal(), abs(z))
end

make_ecdf(Y) = (_ecdf = ecdf(Y); f(x) = _ecdf(x))

function plot_ecdf_pval(F_pval)
    plot(F_pval, 0, 1; label="")
    plot!(identity, 0, 1; label="", ls=:dot, c=:black, alpha=0.5)
    plot!(xtick=0:0.1:1, ytick=0:0.1:1)
    plot!(xguide="nominal significance level α",
        yguide="probability of P-value ≤ α")
    plot!(size=(400, 400))
end

@show pvalue_normal_approx(Binomial(100, 0.5), 60);

@show nulldist = Binomial(100, 1/3)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
@show nulldist = Binomial(100, 1/3)
@show altdist = Binomial(100, 1/2)
X = rand(altdist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
plot(p -> pvalue_normal_approx(Binomial(100, p), 30), 0, 1; label="data: n=100, x=30")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="parameter p", yguide="P-value")

# %%
n, p₀ = 100, 0.4
@show nulldist = Binomial(n, p₀)

p_value_of_null = pvalue_normal_approx(nulldist, 30)

confint95 = find_zeros(p -> pvalue_normal_approx(Binomial(n, p), 30) - 0.05, 0, 1)

o = optimize(p -> -pvalue_normal_approx(Binomial(n, p), 30), 0, 1)
point_estimate = o.minimizer

@show p_value_of_null confint95 point_estimate

plot(p -> pvalue_normal_approx(Binomial(100, p), 30), 0, 1;
    label="P-value function for data: n=100, x=30")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="parameter p", yguide="P-value")
plot!(confint95, fill(0.05, 2); label="95% confidence interval")
scatter!([point_estimate], [1]; label="point estimate")

# %%
@show nulldist = Poisson(30)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
@show nulldist = Poisson(30)
@show altdist = Poisson(45)
X = rand(altdist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
@show nulldist = Poisson(30)

p_value_of_null = pvalue_normal_approx(nulldist, 40)

confint95 = find_zeros(0.01, 100) do λ
    pvalue_normal_approx(Poisson(λ), 40) - 0.05
end

o = optimize(0.01, 100) do λ
    -pvalue_normal_approx(Poisson(λ), 40)
end
point_estimate = o.minimizer

@show p_value_of_null confint95 point_estimate

plot(λ -> pvalue_normal_approx(Poisson(λ), 40), 0.01, 100;
    label="P-value function for data x=40")
plot!(xtick=0:10:100, ytick=0:0.1:1)
plot!(xguide="parameter λ", yguide="P-value")
plot!(confint95, fill(0.05, 2); label="95% confidence interval")
scatter!([point_estimate], [1]; label="point estimate")

# %%
@show nulldist = NegativeBinomial(30, 0.7)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
@show nulldist = NegativeBinomial(30, 0.7)
@show altdist = NegativeBinomial(30, 0.5)
X = rand(altdist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
r, p₀ = 30, 0.7
@show nulldist = NegativeBinomial(r, p₀)

x = 20

p_value_of_null = pvalue_normal_approx(nulldist, x)

confint95 = find_zeros(0.01, 1) do p
    pvalue_normal_approx(NegativeBinomial(r, p), x) - 0.05
end

o = optimize(0.01, 1) do p
    -pvalue_normal_approx(NegativeBinomial(r, p), x)
end
point_estimate = o.minimizer

@show p_value_of_null confint95 point_estimate

plot(p -> pvalue_normal_approx(NegativeBinomial(r, p), x), 0.01, 1;
    label="P-value function for data: r=$r, x=$x")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="parameter p", yguide="P-value")
plot!(confint95, fill(0.05, 2); label="95% confidence interval")
scatter!([point_estimate], [1]; label="point estimate")
plot!(legend=:topleft)

# %%
A = [
    115 90
    80 110
]
a, b, c, d = A'

@show A
@show a, b, c, d
@show (a/b)/(c/d)
println()

@show nulldist = Hypergeometric(a+b, c+d, a+c)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
A = [
    115 90
    80 110
]
a, b, c, d = A'

@show A
@show a, b, c, d
@show (a/b)/(c/d)
println()

@show nulldist = Hypergeometric(a+b, c+d, a+c)
@show altdist = FisherNoncentralHypergeometric(a+b, c+d, a+c, 1.8)
X = rand(altdist, 10^5)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
A = [
    115 90
    80 110
]
a, b, c, d = A'

@show A
@show a, b, c, d
@show (a/b)/(c/d)
println()

@show nulldist = Hypergeometric(a+b, c+d, a+c)
p_value_of_null = pvalue_normal_approx(Hypergeometric(a+b, c+d, a+c), a)

confint95 = find_zeros(0.1, 10) do OR
    pvalue_normal_approx(FisherNoncentralHypergeometric(a+b, c+d, a+c, OR), a) - 0.05
end

o = optimize(0.1, 10) do OR
    -pvalue_normal_approx(FisherNoncentralHypergeometric(a+b, c+d, a+c, OR), a)
end
point_estimate = o.minimizer

@show p_value_of_null confint95 point_estimate

plot(OR -> pvalue_normal_approx(FisherNoncentralHypergeometric(a+b, c+d, a+c, OR), 115), 0.8, 4;
    label="P-value function for data $A")
plot!(xtick=0:0.2:5, ytick=0:0.1:1)
plot!(xguide="parameter OR = odds ratio", yguide="P-value")
plot!(confint95, fill(0.05, 2); label="95% confidence interval")
scatter!([point_estimate], [1]; label="point estimate")

# %%
