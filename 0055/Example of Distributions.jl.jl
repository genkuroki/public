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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# Colab: https://colab.research.google.com/github/genkuroki/public/blob/main/0055/Example%20of%20Distributions.jl.ipynb

# %%
using Pkg: Pkg
haskey(ENV, "COLAB_GPU") && Pkg.add("Distributions")

using Random: rand!
using Distributions
using Plots
default(fmt=:png, legend=false, size=(400, 250))

haskey(ENV, "COLAB_GPU") && run(`apt-get -y install fonts-ipafont-gothic`)
default(fontfamily="ipagp")

"""
    sample_means_naive(dist, n; niters=10^4)

確率分布 `dist` のサイズ `n` のランダムな標本の平均を `niters` 個生成して返す。
ただし、この実装では無駄なメモリアロケーションが発生する。

Example:
```jldoctest
julia> dist = Exponential(1); # 期待値1の指数分布

julia> n = 10; # 標本サイズ10

julia> sample_means_naive(dist, n; niters=10) # 標本平均が10個
10-element Vector{Float64}:
 0.7477545228347106
 1.4361102626055058
 0.9055294429107501
 1.0630035197614673
 0.9325838802514189
 0.6287737443205483
 1.5488765847122619
 1.4094079716586008
 0.3952622382268107
 0.9499173797912887
```
"""
function sample_means_naive(dist, n; niters=10^4)
    Xbar = zeros(niters)
    for i in 1:niters
        X = rand(dist, n) # 無駄なメモリアロケーションが発生
        Xbar[i] = mean(X)
    end
    Xbar
end

"""
    sample_means_single_thread(dist, n; niters=10^4)

確率分布 `dist` のサイズ `n` のランダムな標本の平均を `niters` 個生成して返す。
シングルスレッド版。
"""
function sample_means_single_thread(dist, n; niters=10^4)
    Xbar = zeros(niters)
    Xtmp = zeros(eltype(dist), n)
    for i in 1:niters
        X = rand!(dist, Xtmp) # 無駄なメモリアロケーションの抑制
        Xbar[i] = mean(X)
    end
    Xbar
end

"""
    sample_means(dist, n; niters=10^4)

確率分布 `dist` のサイズ `n` のランダムな標本の平均を `niters` 個生成して返す。
マルチスレッド版。
"""
function sample_means(dist, n; niters=10^4)
    Xbar = zeros(niters)
    nth = Threads.nthreads(:interactive) + Threads.nthreads(:default)
    Xtmp = [zeros(eltype(dist), n) for _ in 1:nth]
    Threads.@threads :static for i in 1:niters
        tid = Threads.threadid()
        X = rand!(dist, Xtmp[tid]) # 無駄なメモリアロケーションの抑制
        Xbar[i] = mean(X)
    end
    Xbar
end

"""
    plot_sample_means(dist, n; niters=10^5)

確率分布 `dist` のサイズ `n` のランダムな標本の平均を `niters` 個生成して、
そのヒストグラムとそれを近似する正規分布を同時プロットする。
"""
function plot_sample_means(dist, n; niters=10^5)
    mu, sigma = mean(dist), std(dist)
    se = sigma / sqrt(n)
    Xbar = sample_means(dist, n; niters)
    histogram(Xbar; norm=true, lc=:match, alpha=0.3)
    plot!(x -> pdf(Normal(mu, se), x), mu - 5se, mu + 5se; lw=1.5)
    plot!(xguide="標本平均")
end

(@doc sample_means_naive) |> display; println()
(@doc sample_means_single_thread) |> display
(@doc sample_means) |> display
@doc plot_sample_means

# %%
# sample_means_naiveでは無駄なメモリアロケーションが発生して遅くなる。
@time sample_means_naive(Exponential(1), 10; niters=10^7);
@time sample_means_naive(Exponential(1), 10; niters=10^7);
@time sample_means_naive(Exponential(1), 10; niters=10^7);
@time sample_means_naive(Exponential(1), 10; niters=10^7);
@time sample_means_naive(Exponential(1), 10; niters=10^7);

# %%
# シングルスレッド版
@time sample_means_single_thread(Exponential(1), 10; niters=10^7);
@time sample_means_single_thread(Exponential(1), 10; niters=10^7);
@time sample_means_single_thread(Exponential(1), 10; niters=10^7);
@time sample_means_single_thread(Exponential(1), 10; niters=10^7);
@time sample_means_single_thread(Exponential(1), 10; niters=10^7);

# %%
# マルチスレッド版
@time sample_means(Exponential(1), 10; niters=10^7);
@time sample_means(Exponential(1), 10; niters=10^7);
@time sample_means(Exponential(1), 10; niters=10^7);
@time sample_means(Exponential(1), 10; niters=10^7);
@time sample_means(Exponential(1), 10; niters=10^7);

# %%
Xbar = sample_means_naive(Exponential(1), 10; niters=10) # 期待値1の指数分布のサイズ10の標本平均を10個

# %%
Xbar = sample_means_single_thread(Exponential(1), 10; niters=10) # 期待値1の指数分布のサイズ10の標本平均を10個

# %%
Xbar = sample_means(Exponential(1), 10; niters=10) # 期待値1の指数分布のサイズ10の標本平均を10個

# %%
plot_sample_means(Normal(), 5)

# %%
plot_sample_means(Uniform(), 5)

# %%
plot_sample_means(Exponential(1), 5)

# %%
plot_sample_means(MixtureModel([Normal(), Normal(20)], [0.95, 0.05]), 20)

# %%
