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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %%
ENV["LINES"] = 200
ENV["COLUMNS"] = 200
using Distributions
using Roots
using StatsPlots
default(fmt=:png)

# %%
logoddsrat(p0, p1) = (log(p0) - log(1-p0)) - (log(p1) - log(1-p1))
logcprat(p0, p1) = -(log(1-p0) - log(1-p1))

function loglikrat(k, n, p0::Real, p1::Real)
    loglikelihood(Binomial(n, p0), k) - loglikelihood(Binomial(n, p1), k)
end

loglikrat(k, n, prior0::Dirac, prior1::Dirac) =
    loglikrat(k, n, mean(prior0), mean(prior1))

loglikrat(k, n, prior0::Beta, prior1::Beta) =
    loglikelihood(BetaBinomial(n, params(prior0)...), k) - loglikelihood(BetaBinomial(n, params(prior1)...), k)

lim_k_C(C, n, p0, p1) = (log(C) + n*logcprat(p0, p1)) / logoddsrat(p0, p1)

function alpha_k_C(k_C, n, p0, p1)
    bin0 = Binomial(n, p0)
    p0 < p1 && return ccdf(bin0, floor(Int, k_C))
    p0 > p1 && return cdf(bin0, ceil(Int, k_C)-1)
    NaN
end

alpha_C(C, n, p0, p1) = alpha_k_C(lim_k_C(C, n, p0, p1), n, p0, p1)

r(x) = round(x; sigdigits=4)

function sim_likrat_test_C(C, n, p0, p1; p=p0, niters=10^7)
    logC = log(C)
    bin = Binomial(n, p)
    c = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:niters
        tid = Threads.threadid()
        k = rand(bin)
        @inbounds c[tid] += loglikrat(k, n, p0, p1) < logC
    end
    sum(c)/niters
end

function sim_likrat_test_k_C(k_C, n, p0, p1; p=p0, niters=10^7)
    bin = Binomial(n, p)
    c = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:niters
        tid = Threads.threadid()
        k = rand(bin)
        c[tid] += p0 < p1 ? k > k_C : p0 > p1 ? k < k_C : NaN
    end
    sum(c)/niters
end

function sim_likrat_test_alpha(α, n, p0, p1; p=p0, niters=10^7)
    bin0 = Binomial(n, p0)
    bin = Binomial(n, p)
    c = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:niters
        tid = Threads.threadid()
        k = rand(bin)
        c[tid] += p0 < p1 ? ccdf(bin0, k) < α : p0 > p1 ? cdf(bin0, k-1) < α : NaN
    end
    sum(c)/niters
end

# %%
n = 100
p0 = 0.3
p1 = 0.5

alphas = [(k_C, alpha_k_C(k_C, n, p0, p1)) for k_C in 35:40]
@show collect(zip(first.(alphas), r.(last.(alphas))))

println()

C = 2.0
@show C
@show alpha_C(C, n, p0, p1) |> r
@show sim_likrat_test_C(C, n, p0, p1)

println()

@show k_C = lim_k_C(C, n, p0, p1) |> r
@show alpha_k_C(k_C, n, p0, p1) |> r
@show sim_likrat_test_k_C(k_C, n, p0, p1)

println()

@show α = 0.0339
@show sim_likrat_test_alpha(α, n, p0, p1)
@show sim_likrat_test_alpha(α, n, p0, p1; p=0.4)
@show sim_likrat_test_alpha(α, n, p0, p1; p=0.5)
;

# %%
n = 100
p0 = 0.3
p1 = 0.1

alphas = [(k_C, alpha_k_C(k_C, n, p0, p1)) for k_C in 20:25]
@show collect(zip(first.(alphas), r.(last.(alphas))))

println()

C = 100.0
@show C
@show alpha_C(C, n, p0, p1) |> r
@show sim_likrat_test_C(C, n, p0, p1)

println()

@show k_C = lim_k_C(C, n, p0, p1) |> r
@show alpha_k_C(k_C, n, p0, p1) |> r
@show sim_likrat_test_k_C(k_C, n, p0, p1)

println()

@show α = 0.0478
@show sim_likrat_test_alpha(α, n, p0, p1)
@show sim_likrat_test_alpha(α, n, p0, p1; p=0.2)
@show sim_likrat_test_alpha(α, n, p0, p1; p=0.1)
;

# %%
n = 100
p0 = 0.3

alphas = [(C, p1, alpha_C(C, n, p0, p1)) for C in (0.1, 0.33, 1.0, 3.0, 10.0, 30.0, 100.0), p1 in (0.31, 0.4, 0.5, 0.6, 0.7)]

# %%
n = 100
p0 = 0.3

alphas = [(C, p1, alpha_C(C, n, p0, p1)) for C in (0.1, 0.33, 1.0, 3.0, 10.0, 30.0, 100.0), p1 in (0.25, 0.20, 0.15, 0.10, 0.05)]

# %%
p0 = 0.3
C = 2.0

@show p0, C
alphas = [(n, p1, alpha_C(C, n, p0, p1)) for n in (10, 30, 100, 300, 1000, 3000, 10000), p1 in (0.301, 0.4, 0.5, 0.6, 0.7)]

# %%
p0 = 0.3
C = 10.0

@show p0, C
alphas = [(n, p1, alpha_C(C, n, p0, p1)) for n in (10, 30, 100, 300, 1000, 3000, 10000), p1 in (0.301, 0.4, 0.5, 0.6, 0.7)]

# %%
chareq(a, b, p, t) = t==1 ? p*(a+b)-b : (p*t^a - 1 + (1-p)*t^(-b)) / (t-1)

function sol_chareq(a, b, p; alg=A42())
    f(t) = chareq(a, b, p, t)
    if p < b/(a+b)
        find_zero(f, (1.0, 1e6), alg)
    elseif p > b/(a+b)
        find_zero(f, (1e-6, 1.0), alg)
    else
        one(p)
    end
end

probapprox(p, q) = q ≤ 1/2 ? p ≈ q : 1-p ≈ 1-q

function prob_ge_M(a, b, p, M, N, x=0; alg=A42())
    probapprox(p, b/(a+b)) && return N/(M+N)
    λ = sol_chareq(a, b, p; alg)
    (λ^(x+N) - 1) / (λ^(M+N) - 1)
end

function random_walk(a, b, p, M, N, x=0.0; nmax=10^3)
    X = Float64[]
    (a ≤ 0 || b ≤ 0) && return X
    for n in 1:nmax
        x += rand() < p ? a : -b
        push!(X, x)
        (x ≤ -N || M ≤ x) && return X
    end
    X
end

function sim_random_walks(a, b, p, M, N, x=0.0; nmax=10^3, niters=10^4)
    X = Vector{Vector{Float64}}(undef, niters)
    for i in 1:niters
        X[i] = random_walk(a, b, p, M, N, x; nmax)
    end
    X, last.(X)
end

function _abp(prior0, prior1, prior_true)
    p0, p1, p_true = mean(prior0), mean(prior1), mean(prior_true)
    logOR, logA = logoddsrat(p0, p1), logcprat(p0, p1)
    if p0 == p1
        a, b, p = 0.0, 0.0, p_true
    elseif p0 < p1 # logOR < logA < 0
        a, b, p = -logA, logA-logOR, 1-p_true
    else  # 0 < logA < logOR
        a, b, p = logOR-logA, logA, p_true
    end
    a, b, p
end

_MN(C0, C1) = log(C1), -log(C0)

function _abpMN(p0, p1, p_true, C0, C1)
    a, b, p = _abp(p0, p1, p_true)
    M, N = _MN(C0, C1)
    a, b, p, M, N
end

function prob_accept(p0, p1, p_true, C0, C1)
    p0 == p1 && return NaN
    a, b, p, M, N = _abpMN(p0, p1, p_true, C0, C1)
    prob_ge_M(a, b, p, M, N)
end

_rand_p_true(p::Real) = p
_rand_p_true(prior::UnivariateDistribution) = rand(prior)

function random_likrat_test(p0, p1, prior_true, C0, C1; nmax=10^3)
    logC0, logC1 = log(C0), log(C1)
    LLR = Float64[]
    k =0
    for n in 1:nmax
        p_true = _rand_p_true(prior_true)
        k += rand(Bernoulli(p_true))
        llr = loglikrat(k, n, p0, p1)
        push!(LLR, llr)
        (llr ≤ logC0 || logC1 ≤ llr) && return LLR
    end
    LLR
end

function sim_likrat_tests(p0, p1, p_true, C0, C1; nmax=10^3, niters=10^4)
    LLR = Vector{Vector{Float64}}(undef, niters)
    for i in 1:niters
        LLR[i] = random_likrat_test(p0, p1, p_true, C0, C1; nmax)
    end
    LLR, last.(LLR), length.(LLR)
end

# %%
a, b, p = 3, 2, 0.6
plot(t -> chareq(a, b, p, t), 0.6, 1)

# %%
a, b, p, M, N = 0.4, 0.5, 0.55, 5, 5
@show a b p M N
@show prob_ge_M(a, b, p, M, N)
@time X, lastX = sim_random_walks(a, b, p, M, N; nmax=10^4, niters=10^5)
@show alpha1 = mean(lastX .≥ M)
@show alpha0 = mean(lastX .≤ -N)
@show alpha1 + alpha0

plot()
for i in 1:100
    plot!(X[i]; label="", c=1, lw=0.5, alpha=0.5)
end
plot!()

# %%
p0, p1, p, C0, C1 = 0.3, 0.5, 0.3, 1/10, 10
logC0, logC1 = log(C0), log(C1)
@show p0 p1 p C0 C1 logC0 logC1
@show prob_accept(p0, p1, p, C0, C1)
@time LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p, C0, C1; nmax=10^3, niters=10^4)
@show mean(lenLLR) std(lenLLR)
@show alpha1 = mean(lastLLR .≥ logC1)
@show alpha0 = mean(lastLLR .≤ logC0)
@show alpha1 + alpha0

plot()
for i in 1:100
    plot!(LLR[i]; label="", c=1, lw=0.5, alpha=0.5)
end
plot!()

# %%
p0, p1, p, C0, C1 = 0.3, 0.4, 0.3, 1/16, 16
logC0, logC1 = log(C0), log(C1)
@show p0 p1 p C0 C1 logC0 logC1
@show prob_accept(p0, p1, p, C0, C1)
@time LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p, C0, C1; nmax=10^3, niters=10^4)
@show mean(lenLLR) std(lenLLR)
@show alpha1 = mean(lastLLR .≥ logC1)
@show alpha0 = mean(lastLLR .≤ logC0)
@show alpha1 + alpha0

plot()
for i in 1:100
    plot!(LLR[i]; label="", c=1, lw=0.5, alpha=0.5)
end
plot!()

# %%
p0, p1, p_true, C0, C1 = 0.3, 0.5, 0.3, 1/10, 10
logC0, logC1 = log(C0), log(C1)
@show p0, p1, p_true, C0, C1
@show logC0, logC1
@show prob_accept(p0, p1, p_true, C0, C1)
@time LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p, C0, C1; nmax=10^3, niters=10^4)
@show mean(lenLLR) std(lenLLR)
@show alpha1 = mean(lastLLR .≥ logC1)
@show alpha0 = mean(lastLLR .≤ logC0)
@show alpha1 + alpha0

plot()
for i in 1:100
    plot!(LLR[i]; label="", c=1, lw=0.5, alpha=0.5)
end
plot!()

# %% tags=[]
p0, p1, p_true, C0, C1 = 0.3, 0.4, 0.3, 1/20, 10
logC0, logC1 = log(C0), log(C1)
@show p0, p1, p_true, C0, C1
@show logC0, logC1
@show prob_accept(p0, p1, p_true, C0, C1)
@time LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p, C0, C1; nmax=10^3, niters=10^4)
@show mean(lenLLR) std(lenLLR)
@show alpha1 = mean(lastLLR .≥ logC1)
@show alpha0 = mean(lastLLR .≤ logC0)
@show alpha1 + alpha0

plot()
for i in 1:100
    plot!(LLR[i]; label="", c=1, lw=0.5, alpha=0.5)
end
plot!()

# %%
p0, p1, p_true, C0, C1 = 0.5, 0.3, 0.3, 1/10, 16
logC0, logC1 = log(C0), log(C1)
@show p0, p1, p_true, C0, C1
@show logC0, logC1
@show prob_accept(p0, p1, p_true, C0, C1)
@time LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p, C0, C1; nmax=10^3, niters=10^4)
@show mean(lenLLR) std(lenLLR)
@show alpha1 = mean(lastLLR .≥ logC1)
@show alpha0 = mean(lastLLR .≤ logC0)
@show alpha1 + alpha0

plot()
for i in 1:100
    plot!(LLR[i]; label="", c=1, lw=0.5, alpha=0.5)
end
plot!()

# %%
p0, p1, p_true = 0.3, 0.6, 0.6
p = p0 < p1 ? 1-p_true : p_true
a = p0 < p1 ? log(1-p0)-log(1-p1) : log(p0)-log(p1)
b = p0 < p1 ? log(p1)-log(p0) : log(1-p1)-log(1-p0)
@show exp(a), exp(-b), p
λ = sol_chareq(a, b, p)
@show λ, 1/λ

@show p*λ^a + (1-p)*λ^(-b)
;

# %%
p0, p1, p_true = 0.7, 0.3, 0.3
p = p0 < p1 ? 1-p_true : p_true
a = p0 < p1 ? log(1-p0)-log(1-p1) : log(p0)-log(p1)
b = p0 < p1 ? log(p1)-log(p0) : log(1-p1)-log(1-p0)
λ = sol_chareq(a, b, p)
@show λ, 1/λ
;

# %%
@show α, β = 0.05, 0.20
@show M, N = log((1-α)/(1-β)), log(β/α)
@show (exp(-N)-1)/(exp(-M-N)-1), (exp(N)-1)/(exp(M+N)-1)

p0, p1, p_true, C0, C1 = 0.3, 0.35, 0.3, α/β, (1-α)/(1-β)
logC0, logC1 = log(C0), log(C1)
@show p0, p1, p_true, C0, C1
@show logC0, logC1
@show prob_accept(p0, p1, p_true, C0, C1)
@time LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p_true, C0, C1; nmax=10^4, niters=10^4)
@show mean(lenLLR) std(lenLLR)
@show alpha1 = mean(lastLLR .≥ logC1)
@show alpha0 = mean(lastLLR .≤ logC0)
@show alpha1 + alpha0

plot()
for i in 1:100
    plot!(LLR[i]; label="", c=1, lw=0.5, alpha=0.5)
end
hline!([logC0, logC1]; c=2, label="")
plot!()

# %%
@show α, β = 0.05, 0.20
@show M, N = log((1-α)/(1-β)), log(β/α)
@show (exp(-N)-1)/(exp(-M-N)-1), (exp(N)-1)/(exp(M+N)-1)

p0, p1, p_true, C0, C1 = 0.3, 0.35, 0.35, α/β, (1-α)/(1-β)
logC0, logC1 = log(C0), log(C1)
@show p0, p1, p_true, C0, C1
@show logC0, logC1
@show prob_accept(p0, p1, p_true, C0, C1)
@time LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p_true, C0, C1; nmax=10^4, niters=10^4)
@show mean(lenLLR) std(lenLLR)
@show alpha1 = mean(lastLLR .≥ logC1)
@show alpha0 = mean(lastLLR .≤ logC0)
@show alpha1 + alpha0

plot()
for i in 1:100
    plot!(LLR[i]; label="", c=1, lw=0.5, alpha=0.5)
end
hline!([logC0, logC1]; c=2, label="")
plot!()

# %% [markdown]
# ## 尤度比が閾値を超えるまでデータの数値を取得し続ける検定法
#
# $
# \newcommand\op{\operatorname}
# $__データの数値:__ $x_1, x_2, x_3, \ldots$ (各$x_i$は$1$または$0$)
#
# $k=x_1+x_2+\cdots+x_n$とおく.
#
# __モデル:__ 成功率$p$のBernoulli試行
#
# $0<p_0<1$, $0<p_1<1$, $p_0\ne p_1$と仮定する.
#
# __帰無仮説:__ $p=p_0$
#
# __対立仮説:__ $p=p_1$
#
# __尤度比__: $n$番目の尤度比を次のように定める:
# $$
# \op{LR}_n = \frac{p_0^k(1-p_0)^{n-k}}{p_1^k(1-p_1)^{n-k}}.
# $$
#
# __閾値の設定:__ $0<\alpha<1$, $0<\beta<1$, $\alpha+\beta<1$と仮定し, 次のようにおく:
# $$
# C_0 = \frac{\alpha}{1-\beta}, \quad
# C_1 = \frac{1-\alpha}{\beta}.
# $$
#
# __判定法:__ $i=1,2,\ldots,n-1$について$C_0\le\op{LR}_i\le C_1$となり, $i=n$で(初めて)そうならなかったとき,
#
# * $\op{LR}_n < C_0$ならば帰無仮説$p=p_0$を棄却する.
# * $\op{LR}_n > C_1$ならば対立仮説$p=p_1$を棄却する.
#
# __αエラー率とβエラー率__
#
# * データの数値が帰無仮説$p=p_0$の下でのモデルの確率分布で生成されているとき, 帰無仮説が棄却される確率は$\alpha$で近似される.
# * データの数値が帰無仮説$p=p_1$の下でのモデルの確率分布で生成されているとき, 対立仮説が棄却される確率は$\beta$で近似される.

# %%
function plot_tests(; α=0.05, β=0.20, p0=0.3, p1=0.5, p_true=0.3, cc=0.0, nmax=1000, niters=10^5)
    println()
    @show α, β, p0, p1, p_true
    a, b, p = _abp(p0, p1, p_true)
    @show a, b, cc
    @show C0 = α/(1-β) * exp(-cc*b)
    @show C1 = (1-α)/β * exp(cc*a)
    println()
    LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p_true, C0, C1; nmax, niters)
    @show mean_lenLLR = mean(lenLLR)
    @show std_lenLLR = std(lenLLR)
    @show median_lenLLR = median(lenLLR)
    println()
    @show samplesize(β; p0=mean(p0), p1=mean(p1), α)
    @show power(ceil(Int, mean_lenLLR); p0=mean(p0), p1=mean(p1), α)
    println()
    @show prob_reject_H0 = mean(lastLLR .< log(C0))
    @show prob_reject_H1 = mean(lastLLR .> log(C1))
    @show prob_reject_H0 + prob_reject_H1

    plot()
    for i in 1:200
        plot!(0:length(LLR[i]), [0.0; LLR[i]]; label="", lw=1, alpha=0.3)
    end
    hline!([log(C0)]; c=:red, label="log(C0)")
    hline!([log(C1)]; c=:blue, label="log(C1)")
    plot!(legendfontsize=12)
    plot!(xguide="n", yguide="log(LRₙ)", guidefontsize=14)
    plot!(size=(800, 300))
    plot!(leftmargin=4Plots.mm, bottommargin=4Plots.mm)
end

function pvalue(k, n, p0, p1)
    if p0 == p1
        NaN
    elseif p0 < p1
        ccdf(Binomial(n, p0), k-1)
    else
        cdf(Binomial(n, p0), k)
    end
end

function hypothesis_test(k, n, p0, p1, α)
    pvalue(k, n, p0, p1) < α
end

function power(n; p0=0.3, p1=0.5, α=0.05)
    bin1 = Binomial(n, p1)
    sum(hypothesis_test(k, n, p0, p1, α) * pdf(bin1, k) for k in support(bin1))
end

function power_mc(n; p0=0.3, p1=0.5, α=0.05, niters=10^5)
    bin1 = Binomial(n, p1)
    c = 0
    for i in 1:niters
        k = rand(bin1)
        c += hypothesis_test(k, n, p0, p1, α)
    end
    c/niters
end

function samplesize(β; p0=0.3, p1=0.5, α=0.05, powerfunc=power)
    n = 1
    while !(1-β ≤ powerfunc(n; p0, p1, α) && 1-β ≤ powerfunc(n+1; p0, p1, α))
        n += 1
    end
    n
end

loglikrat(p, prior0::Distribution, prior1::Distribution) =
    loglikelihood(prior0, p) - loglikelihood(prior1, p)

# %%
plot_tests(; α=0.05, β=0.20, p0=0.3, p1=0.5, p_true=0.3) |> display
plot_tests(; α=0.05, β=0.20, p0=0.3, p1=0.5, p_true=0.5) |> display

# %%
plot_tests(; α=0.025, β=0.20, p0=0.3, p1=0.5, p_true=0.3) |> display
plot_tests(; α=0.025, β=0.20, p0=0.3, p1=0.5, p_true=0.5) |> display

# %%
plot_tests(; α=0.05, β=0.20, p0=0.3, p1=0.2, p_true=0.3) |> display
plot_tests(; α=0.05, β=0.20, p0=0.3, p1=0.2, p_true=0.2) |> display

# %%
plot_tests(; α=0.025, β=0.20, p0=0.3, p1=0.2, p_true=0.3) |> display
plot_tests(; α=0.025, β=0.20, p0=0.3, p1=0.2, p_true=0.2) |> display

# %%
plot_tests(; α=0.025, β=0.20, p0=0.5, p1=0.6, p_true=0.5) |> display
plot_tests(; α=0.025, β=0.20, p0=0.5, p1=0.6, p_true=0.6) |> display

# %%
plot_tests(; α=0.025, β=0.20, p0=0.4, p1=0.6, p_true=0.35) |> display
plot_tests(; α=0.025, β=0.20, p0=0.4, p1=0.6, p_true=0.4) |> display
plot_tests(; α=0.025, β=0.20, p0=0.4, p1=0.6, p_true=0.5) |> display
plot_tests(; α=0.025, β=0.20, p0=0.4, p1=0.6, p_true=0.6) |> display
plot_tests(; α=0.025, β=0.20, p0=0.4, p1=0.6, p_true=0.7) |> display

# %%
m = 1
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
plot(prior0; label="prior0")
plot!(prior1; label="prior1")
plot!(size=(400, 250)) |> display

@show loglikrat(mean(prior0), prior0, prior1)
@show loglikrat(mean(prior1), prior0, prior1)

nmax = 300
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior0), nmax) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior1), nmax) |> display

# %%
m = 10
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
plot(prior0; label="prior0")
plot!(prior1; label="prior1")
plot!(size=(400, 250)) |> display

@show loglikrat(mean(prior0), prior0, prior1)
@show loglikrat(mean(prior1), prior0, prior1)

nmax = 2000
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior0), nmax) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior1), nmax) |> display

# %%
m = 100
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
plot(prior0; label="prior0")
plot!(prior1; label="prior1")
plot!(size=(400, 250)) |> display

@show loglikrat(mean(prior0), prior0, prior1)
@show loglikrat(mean(prior1), prior0, prior1)

nmax = 2000
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior0), nmax) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior1), nmax) |> display

# %%
m = 1000
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
plot(prior0; label="prior0")
plot!(prior1; label="prior1")
plot!(size=(400, 250)) |> display

@show loglikrat(mean(prior0), prior0, prior1)
@show loglikrat(mean(prior1), prior0, prior1)

nmax = 2000
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior0), nmax) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior1), nmax) |> display

# %%
m = 100
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
@show loglikrat((mean(prior0)+mean(prior1))/2, prior0, prior1)
println()
plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior0)-0.05) |> display
plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior0)) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=(mean(prior0)+mean(prior1))/2) |> display
plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior1)) |> display
plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=mean(prior1)+0.10) |> display

# %%
m = 1
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
plot(prior0; label="prior0")
plot!(prior1; label="prior1")
plot!(size=(400, 250)) |> display

@show loglikrat(mean(prior0), prior0, prior1)
@show loglikrat(mean(prior1), prior0, prior1)

nmax = 300
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=prior0, nmax) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=prior1, nmax) |> display

# %%
m = 10
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
plot(prior0; label="prior0")
plot!(prior1; label="prior1")
plot!(size=(400, 250)) |> display

@show loglikrat(mean(prior0), prior0, prior1)
@show loglikrat(mean(prior1), prior0, prior1)

nmax = 1000
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=prior0, nmax) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=prior1, nmax) |> display

# %%
m = 100
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
plot(prior0; label="prior0")
plot!(prior1; label="prior1")
plot!(size=(400, 250)) |> display

@show loglikrat(mean(prior0), prior0, prior1)
@show loglikrat(mean(prior1), prior0, prior1)

nmax = 2000
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=prior0, nmax) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=prior1, nmax) |> display

# %%
m = 1000
prior0 = Beta(4m, 6m)
prior1 = Beta(6m, 4m)
plot(prior0; label="prior0")
plot!(prior1; label="prior1")
plot!(size=(400, 250)) |> display

@show loglikrat(mean(prior0), prior0, prior1)
@show loglikrat(mean(prior1), prior0, prior1)

nmax = 2000
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=prior0, nmax) |> display
@time plot_tests(; α=0.025, β=0.20, p0=prior0, p1=prior1, p_true=prior1, nmax) |> display

# %%
