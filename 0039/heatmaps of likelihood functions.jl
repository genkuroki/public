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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions
using StatsBase: ecdf
using StatsFuns: logsumexp
using StatsPlots
default(fmt=:png,
    titlefontsize=10, plot_titlefontsize=10,
    guidefontsize=9, legendfontsize=8, tickfontsize=6)

# %% [markdown]
# $$
# p(x|a,b)
# = (1-a)\frac{e^{-x^2/2}}{\sqrt{2\pi}} + a\frac{e^{-(x-b)^2/2}}{\sqrt{2\pi}}
# = \frac{e^{-x^2/2}}{\sqrt{2\pi}}
# \left(
# 1 - a + a e^{bx-b^2/2}
# \right)
# $$
#
# のとき
#
# $$
# \log p(x|a,b) 
# = -\frac{x^2}{2} - \log\sqrt{2\pi}
# + \log\left(
# 1 - a + a e^{bx-b^2/2}
# \right)
# $$

# %%
mixnormal(a, b) = MixtureModel([Normal(), Normal(b, 1)], [1-a, a])

function logpdf_mixnormal(a, b, x)
    -x^2/2 - log(√(2π)) + log1p(-a*(1 - exp(b*x - b^2/2)))
end

function pdf_mixnormal(a, b, x)
    exp(logpdf_mixnormal(a, b, x))
end

function loglik_mixnormal(a, b, X)
    sum(logpdf_mixnormal(a, b, x) for x in X)
end

# %%
X = randn(10)
a, b = 0.3, 4
plot(x -> pdf_mixnormal(a, b, x), -5, 10; label="mixnormal(a=$a, b=$b)")
plot!(x -> pdf(mixnormal(a, b), x), -5, 10; label="", ls=:dash)
plot!(; xtick=-10:10)

# %%
function plot_lik_mixnormal(; a = 0.5, b = 1.0, n = 100)
    PP = Vector{Any}(undef, 20)
    as = range(0, 1, 400)
    bs = range(-4, 4, 400)
    Threads.@threads for i in 1:20
        X = rand(mixnormal(a, b), n)
        logz = loglik_mixnormal.(as, bs', Ref(X))
        maxlogz = maximum(logz)
        w = @. exp(logz - maxlogz)
        P = heatmap(bs, as, w; colorbar=false, tickfontsize=5)
        plot!(xtick=-10:10, ytick=-0:0.1:1)
        PP[i] = P
    end
    plot(PP...; layout=(4, 5), size=(1000, 800))
    plot!(plot_title="true a = $a,  true b = $b,  n = $n")
end

# %%
@time plot_lik_mixnormal(; a = 0.5, b = 1.0, n = 100)

# %%
@time plot_lik_mixnormal(; a = 0.5, b = 0.5, n = 100)

# %%
@time plot_lik_mixnormal(; a = 0.5, b = 0.2, n = 100)

# %%
@time plot_lik_mixnormal(; a = 0.5, b = 0.1, n = 100)

# %%
@time plot_lik_mixnormal(; a = 0.5, b = 0.0, n = 100)

# %%
regtanh(a, b, x) = Normal(b*tanh(a*x), 1)
loglik_regtanh(a, b, X, Y) = sum(logpdf(regtanh(a, b, x), y) for (x, y) in zip(X, Y))

function plog_lik_regtanh(; a=1, b=1, n=100)
    PP = Vector{Any}(undef, 20)
    X = range(-1, 1, n)
    as = range(-3, 3, 400)
    bs = range(-3, 3, 400)
    Threads.@threads for i in 1:20
        Y = [rand(regtanh(a, b, x)) for x in X]
        logz = loglik_regtanh.(as, bs', Ref(X), Ref(Y))
        maxlogz = maximum(logz)
        w = @. exp(logz - maxlogz)
        P = heatmap(bs, as, w; colorbar=false, tickfontsize=5)
        plot!(xtick=-10:10, ytick=-10:10)
        PP[i] = P
    end
    plot(PP...; layout=(4, 5), size=(1000, 800))
    plot!(plot_title="true a = $a,  true b = $b,  n = $n")
end

# %%
@time plog_lik_regtanh(; a=1, b=1, n=100)

# %%
@time plog_lik_regtanh(; a=0.5, b=1, n=100)

# %%
@time plog_lik_regtanh(; a=0.2, b=1, n=100)

# %%
@time plog_lik_regtanh(; a=0.1, b=1, n=100)

# %%
@time plog_lik_regtanh(; a=0, b=1, n=100)

# %%
@time plog_lik_regtanh(; a=1, b=1, n=1000)

# %%
@time plog_lik_regtanh(; a=0.5, b=1, n=1000)

# %%
@time plog_lik_regtanh(; a=0.2, b=1, n=1000)

# %%
@time plog_lik_regtanh(; a=0.1, b=1, n=1000)

# %%
@time plog_lik_regtanh(; a=0, b=1, n=1000)

# %%
