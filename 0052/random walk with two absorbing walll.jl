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
loglikrat(k, n, p0, p1) = loglikelihood(Binomial(n, p0), k) - loglikelihood(Binomial(n, p1), k)
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

function _abpMN(p0, p1, p_true, C0, C1)
    logOR, logA = logoddsrat(p0, p1), logcprat(p0, p1)
    M, N = log(C1), -log(C0)
    if p0 == p1
        a, b, p = 0.0, 0.0, p_true
    elseif p0 < p1 # logOR < logA < 0
        a, b, p = -logA, logA-logOR, 1-p_true
    else  # 0 < logA < logOR
        a, b, p = logOR-logA, logA, p_true
    end
    a, b, p, M, N
end

function prob_accept(p0, p1, p_true, C0, C1)
    p0 == p1 && return NaN
    a, b, p, M, N = _abpMN(p0, p1, p_true, C0, C1)
    prob_ge_M(a, b, p, M, N)
end

function random_likrat_test(p0, p1, p_true, C0, C1; nmax=10^3)
    logC0, logC1 = log(C0), log(C1)
    LLR = Float64[]
    k =0
    for n in 1:nmax
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
#??find_zero

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

# %%
@show α, β = 0.025, 0.20
@show M, N = log((1-α)/(1-β)), log(β/α)
@show (exp(-N)-1)/(exp(-M-N)-1), (exp(N)-1)/(exp(M+N)-1)

println()

p0, p1, p_true, C0, C1 = 0.30, 0.35, 0.30, α/β, (1-α)/(1-β)
logC0, logC1 = log(C0), log(C1)
@show p0, p1, p_true, C0, C1
@show logC0, logC1
@show _abpMN(p0, p1, p_true, M, N)
@show prob_accept(p0, p1, p_true, C0, C1)

println()

@time LLR, lastLLR, lenLLR = sim_likrat_tests(p0, p1, p_true, C0, C1; nmax=10^4, niters=10^4)
@show mean(lenLLR) std(lenLLR)
@show alpha1 = mean(lastLLR .≥ logC1)
@show alpha0 = mean(lastLLR .≤ logC0)
@show alpha1 + alpha0

plot()
for i in 1:200
    plot!(LLR[i]; label="", c=1, lw=1, alpha=0.2)
end
hline!([logC0, logC1]; c=2, label="")
plot!()

# %%
