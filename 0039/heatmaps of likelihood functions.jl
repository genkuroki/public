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
# \begin{aligned}
# p(x|a,b)
# &= (1-a)\frac{e^{-x^2/2}}{\sqrt{2\pi}} + a\frac{e^{-(x-b)^2/2}}{\sqrt{2\pi}}
# \\ &= \frac{e^{-x^2/2}}{\sqrt{2\pi}}
# \left(
# 1 - a + a e^{bx-b^2/2}
# \right)
# \\ &= \frac{e^{-x^2/2}}{\sqrt{2\pi}}
# \left(
# 1 + a\left(e^{bx-b^2/2} - 1\right)
# \right)
# \end{aligned}
# $$
#
# のとき
#
# $$
# \log p(x|a,b) 
# = -\frac{x^2}{2} - \log\sqrt{2\pi}
# + \log\left(
# 1 + a\left(e^{bx-b^2/2} - 1\right)
# \right).
# $$
#
# データ $X=(x_1,\ldots,x_n)$ の対数尤度:
#
# $$
# L(a,b|X) = \sum_{i=1}^n \log p(x_i|a,b).
# $$

# %%
?log1p

# %%
?expm1

# %%
mixnormal(a, b) = MixtureModel([Normal(), Normal(b, 1)], [1-a, a])

function logpdf_mixnormal(a, b, x)
    -x^2/2 - log(√(2π)) + log1p(a*expm1(b*x - b^2/2))
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
    as = range(0, 1, 400)
    bs = range(-4, 4, 400)
    ws = Vector{Matrix{Float64}}(undef, 20)
    Threads.@threads for i in 1:20
        X = rand(mixnormal(a, b), n)
        logz = loglik_mixnormal.(as, bs', Ref(X))
        maxlogz = maximum(logz)
        ws[i] = @. exp(logz - maxlogz)
    end
    PP = Vector{Any}(undef, 20)
    for i in 1:20
        P = heatmap(bs, as, ws[i]; colorbar=false, tickfontsize=5)
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

# %% [markdown]
# モデル:
#
# $$
# p(y|a,b,x) = \frac{1}{\sqrt{2\pi}}\exp\left(\frac{(y - b\tanh ax)^2}{2}\right).
# $$
#
# データ $X=(x_1,\ldots,x_n)$, $Y=(y_1,\ldots,y_n)$ の対数尤度:
#
# $$
# \log L(a,b|X,Y) = \sum_{i=1}^n \log p(y_i|a,b,x_i).
# $$

# %%
regtanh(a, b, x) = Normal(b*tanh(a*x), 1)
loglik_regtanh(a, b, X, Y) = sum(logpdf(regtanh(a, b, x), y) for (x, y) in zip(X, Y))

# %%
a, b, n = 3, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)

# %%
a, b, n = 1, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)

# %%
a, b, n = 0.5, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)

# %%
a, b, n = 0.2, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)

# %%
a, b, n = 0.1, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)

# %%
a, b, n = 0, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)

# %%
function plog_lik_regtanh(; a=1, b=1, n=100)
    X = range(-1, 1, n)
    as = range(-3, 3, 400)
    bs = range(-3, 3, 400)
    ws = Vector{Matrix{Float64}}(undef, 20)
    Threads.@threads for i in 1:20
        Y = [rand(regtanh(a, b, x)) for x in X]
        logz = loglik_regtanh.(as, bs', Ref(X), Ref(Y))
        maxlogz = maximum(logz)
        ws[i] = @. exp(logz - maxlogz)
    end
    PP = Vector{Any}(undef, 20)
    for i in 1:20
        P = heatmap(bs, as, ws[i]; colorbar=false, tickfontsize=5)
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
