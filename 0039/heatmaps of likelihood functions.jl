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

# %% tags=[]
?log1p

# %% tags=[]
?expm1

# %% [markdown]
# モデル:
#
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
# 1 + a\left(e^{b(x-b/2)} - 1\right)
# \right).
# \end{aligned}
# $$
#
# このとき
#
# $$
# \log p(x|a,b) 
# = -\frac{x^2}{2} - \log\sqrt{2\pi}
# + \log\left(
# 1 + a\left(e^{b(x-b/2)} - 1\right)
# \right).
# $$
#
# データ $X=(x_1,\ldots,x_n)$ の対数尤度:
#
# $$
# L(a,b|X) = \sum_{i=1}^n \log p(x_i|a,b).
# $$

# %%
mixnormal(a, b) = MixtureModel([Normal(), Normal(b, 1)], [1-a, a])

function logpdf_mixnormal(a, b, x)
    -x^2/2 - log(√(2π)) + log1p(a*expm1(b*(x - b/2)))
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

# %%
function anim_lik_mixnormal(; a = 0.5, b = 0.2, n = 100, nrepl=20, nframes=400)
    as = range(0, 1, 400)
    bs = range(-4, 4, 400)
    Xall = rand(mixnormal(a, b), n + nframes*nrepl)
    ws = Vector{Matrix{Float64}}(undef, nframes+1)
    Threads.@threads for i in 0:nframes
        X = @view Xall[i*nrepl+1:i*nrepl+n]
        logz = loglik_mixnormal.(as, bs', Ref(X))
        maxlogz = maximum(logz)
        ws[1+i] = @. exp(logz - maxlogz)
    end
    @gif for i in [fill(0, 20); 0:nframes; fill(nframes, 20)]
        P = heatmap(bs, as, ws[1+i]; colorbar=false, tickfontsize=5)
        plot!(xtick=-10:10, ytick=-0:0.1:1)
        plot!(xguide="b", yguide="a")
        plot!(plot_title="true a = $a,  true b = $b,  n = $n")
        plot!(size=(400, 400))
    end
end

# %%
@time anim_lik_mixnormal(; a = 0.5, b = 0.2, n = 100)

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
title!("a = $a,  b = $b,  n = $n")

# %%
a, b, n = 1, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

# %%
a, b, n = 0.5, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

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
title!("a = $a,  b = $b,  n = $n")

# %%
a, b, n = 0, 1, 100
X = range(-1, 1, n)
Y = [rand(regtanh(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> b*tanh(a*x), -1, 1; label="$(b==1 ? "" : "$b ")tanh($(a)x)")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

# %%
function plot_lik_regtanh(; a=1, b=1, n=100)
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
@time plot_lik_regtanh(; a=1, b=1, n=100)

# %%
@time plot_lik_regtanh(; a=0.5, b=1, n=100)

# %%
@time plot_lik_regtanh(; a=0.2, b=1, n=100)

# %%
@time plot_lik_regtanh(; a=0.1, b=1, n=100)

# %%
@time plot_lik_regtanh(; a=0, b=1, n=100)

# %%
@time plot_lik_regtanh(; a=1, b=1, n=1000)

# %%
@time plot_lik_regtanh(; a=0.5, b=1, n=1000)

# %%
@time plot_lik_regtanh(; a=0.2, b=1, n=1000)

# %%
@time plot_lik_regtanh(; a=0.1, b=1, n=1000)

# %%
@time plot_lik_regtanh(; a=0, b=1, n=1000)

# %%
function anim_lik_regtanh(; a = 0.2, b = 1.0, n = 100, nrepl=20, nframes=400)
    X = range(-1, 1, n)
    as = range(-3, 3, 400)
    bs = range(-3, 3, 400)
    Yall = Vector{Vector{Float64}}(undef, nframes+1)
    Y = [rand(regtanh(a, b, x)) for x in X]
    Yall[1] = copy(Y)
    for i in 1:nframes
        n0 = mod((i-1)*nrepl, n)
        Y[n0+1:n0+nrepl] .= (rand(regtanh(a, b, x)) for x in @view X[n0+1:n0+nrepl])
        Yall[i+1] = copy(Y)
    end
    ws = Vector{Matrix{Float64}}(undef, nframes+1)
    Threads.@threads for i in 0:nframes
        Y = Yall[i+1]
        logz = loglik_regtanh.(as, bs', Ref(X), Ref(Y))
        maxlogz = maximum(logz)
        ws[i+1] = @. exp(logz - maxlogz)
    end
    @gif for i in [fill(0, 20); 0:nframes; fill(nframes, 20)]
        P = heatmap(bs, as, ws[i+1]; colorbar=false, tickfontsize=5)
        plot!(xtick=-10:10, ytick=-10:10)
        plot!(xguide="b", yguide="a")
        plot!(plot_title="true a = $a,  true b = $b,  n = $n")
        plot!(size=(400, 400))
    end
end

# %%
@time anim_lik_regtanh(; a=0.2, b=1, n=100)

# %% [markdown]
# モデル:
#
# $$
# p(y|a,b,x) = \frac{1}{\sqrt{2\pi}}\exp\left(\frac{1}{2}\left(y - \mu(a,b,x)\right)^2\right),
# \quad
# \mu(a,b,x) = \tanh(ax)\tanh((a-2b)x)\tanh((a+2b)x).
# $$

# %%
modeltanh3(a, b, x) = Normal(tanh(a*x)*tanh((a-2b)*x)*tanh((a+2b)*x), 1)
loglik_modeltanh3(a, b, X, Y) = sum(logpdf(modeltanh3(a, b, x), y) for (x, y) in zip(X, Y))

# %%
a, b, n = 2, 2, 100
X = range(-1, 1, n)
Y = [rand(modeltanh3(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> tanh(a*x)*tanh((a-2b)*x)*tanh((a+2b)*x), -1, 1; label="")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

# %%
a, b, n = 1, 1, 100
X = range(-1, 1, n)
Y = [rand(modeltanh3(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> tanh(a*x)*tanh((a-2b)*x)*tanh((a+2b)*x), -1, 1; label="")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

# %%
a, b, n = 0.5, 0.5, 100
X = range(-1, 1, n)
Y = [rand(modeltanh3(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> tanh(a*x)*tanh((a-2b)*x)*tanh((a+2b)*x), -1, 1; label="")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

# %%
a, b, n = 0.2, 0.2, 100
X = range(-1, 1, n)
Y = [rand(modeltanh3(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> tanh(a*x)*tanh((a-2b)*x)*tanh((a+2b)*x), -1, 1; label="")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

# %%
a, b, n = 0.1, 0.1, 100
X = range(-1, 1, n)
Y = [rand(modeltanh3(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> tanh(a*x)*tanh((a-2b)*x)*tanh((a+2b)*x), -1, 1; label="")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

# %%
a, b, n = 0, 0, 100
X = range(-1, 1, n)
Y = [rand(modeltanh3(a, b, x)) for x in X]
plot(X, Y; label="data")
plot!(; xguide="X", yguide="Y")
plot!(x -> tanh(a*x)*tanh((a-2b)*x)*tanh((a+2b)*x), -1, 1; label="")
plot!(legend=:topleft)
title!("a = $a,  b = $b,  n = $n")

# %%
function plot_lik_modeltanh3(; a=1, b=1, n=100)
    X = range(-1, 1, n)
    as = range(-4, 4, 400)
    bs = range(-4, 4, 400)
    ws = Vector{Matrix{Float64}}(undef, 20)
    Threads.@threads for i in 1:20
        Y = [rand(regtanh(a, b, x)) for x in X]
        logz = loglik_modeltanh3.(as, bs', Ref(X), Ref(Y))
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
@time plot_lik_modeltanh3(; a=1, b=1, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.9, b=0.9, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.8, b=0.8, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.7, b=0.7, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.6, b=0.6, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.5, b=0.5, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.4, b=0.4, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.3, b=0.3, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.2, b=0.2, n=100)

# %%
@time plot_lik_modeltanh3(; a=0.1, b=0.1, n=100)

# %%
@time plot_lik_modeltanh3(; a=0, b=0, n=100)

# %%
function anim_lik_modeltanh3(; a = 0.7, b = 0.7, n = 100, nrepl=20, nframes=400)
    X = range(-1, 1, n)
    as = range(-3, 3, 400)
    bs = range(-3, 3, 400)
    Yall = Vector{Vector{Float64}}(undef, nframes+1)
    Y = [rand(modeltanh3(a, b, x)) for x in X]
    Yall[1] = copy(Y)
    for i in 1:nframes
        n0 = mod((i-1)*nrepl, n)
        Y[n0+1:n0+nrepl] .= (rand(modeltanh3(a, b, x)) for x in @view X[n0+1:n0+nrepl])
        Yall[i+1] = copy(Y)
    end
    ws = Vector{Matrix{Float64}}(undef, nframes+1)
    Threads.@threads for i in 0:nframes
        Y = Yall[i+1]
        logz = loglik_modeltanh3.(as, bs', Ref(X), Ref(Y))
        maxlogz = maximum(logz)
        ws[i+1] = @. exp(logz - maxlogz)
    end
    @gif for i in [fill(0, 20); 0:nframes; fill(nframes, 20)]
        P = heatmap(bs, as, ws[i+1]; colorbar=false, tickfontsize=5)
        plot!(xtick=-10:10, ytick=-10:10)
        plot!(xguide="b", yguide="a")
        plot!(plot_title="true a = $a,  true b = $b,  n = $n")
        plot!(size=(400, 400), titlefontsize=8)
    end
end

# %%
@time anim_lik_modeltanh3(; a=0.7, b=0.7, n=100)

# %%
