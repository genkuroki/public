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

# %% [markdown]
# https://x.com/Amatusora_heiwa/status/1673987379504480261

# %%
using Distributions
using StaticArrays
using StatsPlots
default(fmt=:png)

# %%
A = @SMatrix [
     4  7 10
     6  4  6
    13  3  2
]

# %%
degree_of_freedom(A) = (size(A, 1) - 1)*(size(A, 2) - 1)

function pearson_expectation_value(A)
    n = sum(A)
    c = sum(A; dims=1)
    r = sum(A; dims=2)
    r*c/n
end

function pearson_chisq(A)
    E = pearson_expectation_value(A)
    sum((a - e)^2/e for (a, e) in zip(A, E))
end

function pvalue_pearson_chisq(A)
    df = degree_of_freedom(A)
    chi2 = pearson_chisq(A)
    ccdf(Chisq(df), chi2)
end

function pearson_adjusted_standardized_residual(A)
    n = sum(A)
    c = sum(A; dims=1)
    r = sum(A; dims=2)
    p = c/n
    q = r/n
    E = r*c/n
    (A - E) ./ .√(E .* ((1 .- q)*(1 .- p)))
end

@show pearson_chisq(A)
@show degree_of_freedom(A)
@show pvalue_pearson_chisq(A)
asr = pearson_adjusted_standardized_residual(A)
display(asr)
pval_asr = @. 2ccdf(Normal(), abs(asr))

# %%
c = sum(A; dims=1)
r = sum(A; dims=2)
n = sum(A)
E = r * c / n
prodpoi = product_distribution((Poisson(λ) for λ in vec(E))...)
L = 10^6
X = [reshape(rand(prodpoi), 3, 3) for _ in 1:L]
Chi2 = pearson_chisq.(X)
df = (length(c) - 1)*(length(r) - 1)

P = stephist(Chi2; norm=true, label="Pearson's χ²-statistic")
plot!(Chisq(df); label="Chisq($df)")
plot!(xlim=(-1, 24))
plot!(xguide="chi-squared value", yguide="probability density")

Pval = ccdf.(Chisq(df), Chi2)
Q = plot(α -> count(<(α), Pval)/length(Pval), 0, 0.1; label="")
plot!(identity, 0, 0.1; ls=:dot, label="")
plot!(xguide="α", yguide="probability of P-value < α")

plot(P, Q; size=(1000, 400), layout=@layout [a b{0.4w}])
plot!(bottommargin=4Plots.mm, leftmargin=4Plots.mm)

# %%
ASR = pearson_adjusted_standardized_residual.(X)
PP = []
for i in 1:3, j in 1:3
    P = stephist(getindex.(ASR, i, j); norm=true, bin=50, label="")
    plot!(Normal(); label="")
    push!(PP, P)
end
plot(PP...; size=(1000, 600), layout=(3, 3))

# %%
c = sum(A; dims=1)
r = sum(A; dims=2)
n = sum(A)
p = c / n
q = r / n
P = q * p
mult = Multinomial(n, vec(P))
L = 10^6
X = [reshape(rand(mult), 3, 3) for _ in 1:L]
Chi2 = pearson_chisq.(X)
df = (length(c) - 1)*(length(r) - 1)

P = stephist(Chi2; norm=true, label="Pearson's χ²-statistic")
plot!(Chisq(df); label="Chisq($df)")
plot!(xlim=(-1, 24))
plot!(xguide="chi-squared value", yguide="probability density")

Pval = ccdf.(Chisq(df), Chi2)
Q = plot(α -> count(<(α), Pval)/length(Pval), 0, 0.1; label="")
plot!(identity, 0, 0.1; ls=:dot, label="")
plot!(xguide="α", yguide="probability of P-value < α")

plot(P, Q; size=(1000, 400), layout=@layout [a b{0.4w}])
plot!(bottommargin=4Plots.mm, leftmargin=4Plots.mm)

# %%
ASR = pearson_adjusted_standardized_residual.(X)
PP = []
for i in 1:3, j in 1:3
    P = stephist(getindex.(ASR, i, j); norm=true, bin=50, label="")
    plot!(Normal(); label="")
    push!(PP, P)
end
plot(PP...; size=(1000, 600), layout=(3, 3))

# %%
c = sum(A; dims=1)
r = sum(A; dims=2)
n = sum(A)
p = c / n
q = r / n
P = q * p
mult = [Multinomial(m, vec(p)) for m in vec(r)]
L = 10^6
X = [[rand(mult[1])'; rand(mult[2])'; rand(mult[3])'] for _ in 1:L]
Chi2 = pearson_chisq.(X)
df = (length(c) - 1)*(length(r) - 1)

P = stephist(Chi2; norm=true, bin=200, label="Pearson's χ²-statistic")
plot!(Chisq(df); label="Chisq($df)")
plot!(xlim=(-1, 24))
plot!(xguide="chi-squared value", yguide="probability density")

Pval = ccdf.(Chisq(df), Chi2)
Q = plot(α -> count(<(α), Pval)/length(Pval), 0, 0.1; label="")
plot!(identity, 0, 0.1; ls=:dot, label="")
plot!(xguide="α", yguide="probability of P-value < α")

plot(P, Q; size=(1000, 400), layout=@layout [a b{0.4w}])
plot!(bottommargin=4Plots.mm, leftmargin=4Plots.mm)

# %%
ASR = pearson_adjusted_standardized_residual.(X)
PP = []
for i in 1:3, j in 1:3
    P = stephist(getindex.(ASR, i, j); norm=true, bin=50, label="")
    plot!(Normal(); label="")
    push!(PP, P)
end
plot(PP...; size=(1000, 600), layout=(3, 3))

# %%
c = sum(A; dims=1)
r = sum(A; dims=2)
n = sum(A)
p = c / n
q = r / n
P = q * p
mult = [Multinomial(m, vec(q)) for m in vec(c)]
L = 10^6
X = [[rand(mult[1]) rand(mult[2]) rand(mult[3])] for _ in 1:L]
Chi2 = pearson_chisq.(X)
df = (length(c) - 1)*(length(r) - 1)

P = stephist(Chi2; norm=true, bin=200, label="Pearson's χ²-statistic")
plot!(Chisq(df); label="Chisq($df)")
plot!(xlim=(-1, 24))
plot!(xguide="chi-squared value", yguide="probability density")

Pval = ccdf.(Chisq(df), Chi2)
Q = plot(α -> count(<(α), Pval)/length(Pval), 0, 0.1; label="")
plot!(identity, 0, 0.1; ls=:dot, label="")
plot!(xguide="α", yguide="probability of P-value < α")

plot(P, Q; size=(1000, 400), layout=@layout [a b{0.4w}])
plot!(bottommargin=4Plots.mm, leftmargin=4Plots.mm)

# %%
ASR = pearson_adjusted_standardized_residual.(X)
PP = []
for i in 1:3, j in 1:3
    P = stephist(getindex.(ASR, i, j); norm=true, bin=50, label="")
    plot!(Normal(); label="")
    push!(PP, P)
end
plot(PP...; size=(1000, 600), layout=(3, 3))

# %%
