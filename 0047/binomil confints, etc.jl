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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Roots

function pvalue_wald(k, n, p=0.5)
    p̂ = k/n
    z = (p̂ - p)/√(p̂*(1-p̂)/n)
    2ccdf(Normal(), abs(z))
end

function ci_wald(k, n, α=0.05)
    c = cquantile(Normal(), α/2)
    p̂ = k/n
    L = p̂ - c * √(p̂*(1-p̂)/n)
    U = p̂ + c * √(p̂*(1-p̂)/n)
    [L, U]
end

function ci_wald_fake(k, n, α=0.05)
    p̂ = k/n
    bin = Binomial(n, p̂)
    normal = Normal(mean(bin), std(bin))
    l, u = quantile.(normal, (α/2, 1-α/2))
    [l/n, u/n]
end

function pvalue_wilson(k, n, p=0.5)
    p̂ = k/n
    z = (p̂ - p)/√(p*(1-p)/n)
    2ccdf(Normal(), abs(z))
end

function ci_wilson(k, n, α=0.05)
    c = cquantile(Normal(), α/2)
    p̂ = k/n
    L = 1/(1+c^2/n) * (p̂ + c^2/(2n) - c*√(p̂*(1-p̂)/n + c^2/(4n^2)))
    U = 1/(1+c^2/n) * (p̂ + c^2/(2n) + c*√(p̂*(1-p̂)/n + c^2/(4n^2)))
    [L, U]
end

function ci_wilson2(k, n, α=0.05)
    c = cquantile(Normal(), α/2)
    p̂ = k/n
    A, B, C = 1+c^2/n, p̂+c^2/(2n), p̂^2
    # Ap² - 2Bp + C = 0 を解く.
    sqrtD = √(B^2 - A*C)
    L = (B - sqrtD)/A
    U = (B + sqrtD)/A
    [L, U]
end

function ci_wilson3(k, n, α=0.05)
    p̂ = k/n
    c = cquantile(Normal(), α/2)
    f(p) = (p̂-p)^2 - c^2*p*(1-p)/n
    find_zeros(f, (-0.1, 1.1))
end

# %%
n = 2100
k = 0.16n
@show ci_wald(k, n)
@show ci_wald_fake(k, n)
@show ci_wilson(k, n)
@show ci_wilson2(k, n)
@show ci_wilson3(k, n)
;

# %%
function power(n, p1; p0=0.5, pvaluefunc=pvalue_wald, α=0.05, L=10^8)
    bin1 = Binomial(n, p1)
    c = 0
    for i in 1:L
        k = rand(bin1)
        c += pvaluefunc(k, n, p0) < α
    end
    c/L
end

# %%
power(40, 0.7)

# %% tags=[]
function sim(; n=40, p1=0.7, p0=0.5, prevalence=1e-4, 
        α=0.05, pvaluefunc=pvalue_wald, L=10^8)
    bin1, bin0 = Binomial(n, p1), Binomial(n, p0)
    a = b = c = d = 0
    for i in 1:L
        if rand() > prevalence
            k = rand(bin0)
            if pvaluefunc(k, n, p0) < α
                a += 1
            else
                b += 1
            end
        else
            k = rand(bin1)
            if pvaluefunc(k, n, p0) < α
                c += 1
            else
                d += 1
            end
        end
    end
    [a b; c d]
end

# %%
@time a, c, b, d = sim()

# %%
c/(a+c)

# %%
power(1000, 0.5437)

# %%
@time a, c, b, d = sim(; n=1000, p1=0.5437, L=10^8)

# %%
c/(a+c)

# %%
f(; α=0.05, β=0.2, prevalence=1e-4) = prevalence*(1-β) / ((1-prevalence)*α + prevalence*(1-β))
g(; α=0.05, β=0.2, prevalence=1e-4) = (1-prevalence)*(1-α) / ((1-prevalence)*(1-α) + prevalence*β)
f(), g()

# %%
f(; α=0.05, β=0.0, prevalence=1e-4), g(; α=0.05, β=1.0, prevalence=1e-4)

# %%
using Distributions

@show n = 2100
@show p̂ = 0.16
@show bin = Binomial(n, p̂)
@show normal = Normal(mean(bin), std(bin))
@show q = [quantile(normal, 0.025), quantile(normal, 1-0.025)]
@show q / n;

# %%
√(2100*0.16*(1-0.16))

# %%
using RCall
R"""prop.test(2100*0.16, 2100, p=0.14, correct=F)"""

# %%
n = 2100
k = 0.16n
@show pvalue_wilson(k, n, 0.14)
@show ci_wilson(k, n);

# %%
ENV["COLUMNS"] = 1000
using RCall
@rlibrary TOSTER

n = 22
x = 2collect(1:n)
y = x .- 1
println("x and y:")
Base.print_matrix(stdout, [x'; y'])
brunner_munzel(x, y; paired=true)

# %%
ENV["COLUMNS"] = 1000
using RCall
@rlibrary TOSTER

n = 21
x = 2collect(1:n)
y = x .- 1
println("x and y:")
Base.print_matrix(stdout, [x'; y'])
brunner_munzel(x, y; paired=true)

# %%
ENV["COLUMNS"] = 1000
using RCall
@rlibrary TOSTER

n = 21
x = 2collect(1:n)
y = x .- 1
println("x and y:")
Base.print_matrix(stdout, [x'; y'])
brunner_munzel(x, y)#; paired=true)

# %%
ENV["COLUMNS"] = 1000
using RCall
@rlibrary TOSTER

n = 1000
x = 2collect(1:n)
y = x .- 1
x = [x; 3n]
y = [y; 3n+1]
#println("x and y:")
#Base.print_matrix(stdout, [x'; y'])
brunner_munzel(x, y; paired=true)

# %%
mu = 0.5
BM1 = [mean((xᵢ < yⱼ) + (xᵢ == yⱼ)/2 for xᵢ in x) for yⱼ in y]
BM2 = [mean((yⱼ < xᵢ) + (yⱼ == xᵢ)/2 for yⱼ in y) for xᵢ in x]
BM3 = BM1 - BM2
#@show BM1 BM2 BM3
@show pd = mean(BM2)
@show v = var(BM3)
t = √n * (pd - mu) / √v

# %%
m = mean(BM3)
(sum(BM3 .^ 2) - n*m^2)/(n-1)

# %%
var(BM3)

# %%
sum((BM3 .- m) .^ 2)/(n-1)

# %%
@rput x y;
@show R"x"
@show R"y";

# %%
R"""
mu = 0.5
n = length(x)
N = length(c(y, x))
rx = rank(c(y, x))
rx1 = rx[1:n]
rx2 = rx[(n+1):N]
rix1 = rank(y)
rix2 = rank(x)
BM1 = 1/n * (rx1 - rix1)
BM2 = 1/n * (rx2 - rix2)
BM3 = BM1 - BM2
pd = mean(BM2)
"""

# %%
R"""
m = mean(BM3)
v = (sum(BM3 ^ 2) - n * m ^ 2) / (n - 1)
"""

# %%
R"""
v0 = (v == 0)
v[v0] = 1/n
v
"""

# %%
R"""
test_stat = sqrt(n) * (pd - mu) / sqrt(v)
"""

# %%
R"""
sum(BM3 ^ 2)
"""

# %%
sum(BM3 .^ 2)

# %%
rcopy(R"""BM3""") - BM3

# %%
