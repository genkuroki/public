# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# # Ten Lines in Julia
#
# * Author: Gen Kuroki
# * Date: 2021-03-16,17, 2021-05-30, 2021-08-05
# * Copyright (c) 2021 Gen Kuroki
# * License: https://opensource.org/licenses/MIT
# * Twitter: https://twitter.com/genkuroki/status/1371793118203322374

# %% [markdown] toc=true
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#probability-that-the-confidence-interval-contains-the-true-value" data-toc-modified-id="probability-that-the-confidence-interval-contains-the-true-value-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>probability that the confidence interval contains the true value</a></span></li><li><span><a href="#distribution-of-confidence-intervals-of-the-binomial-distribution-model" data-toc-modified-id="distribution-of-confidence-intervals-of-the-binomial-distribution-model-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>distribution of confidence intervals of the binomial distribution model</a></span></li><li><span><a href="#Fisher's-exact-test-version" data-toc-modified-id="Fisher's-exact-test-version-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Fisher's exact test version</a></span></li><li><span><a href="#3-dimensional-random-walk" data-toc-modified-id="3-dimensional-random-walk-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>3-dimensional random walk</a></span></li><li><span><a href="#Ising-model-on-the-2-dimensional-square-lattice" data-toc-modified-id="Ising-model-on-the-2-dimensional-square-lattice-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Ising model on the 2-dimensional square lattice</a></span></li><li><span><a href="#Kuramoto-model" data-toc-modified-id="Kuramoto-model-6"><span class="toc-item-num">6&nbsp;&nbsp;</span>Kuramoto model</a></span></li><li><span><a href="#Julia-set" data-toc-modified-id="Julia-set-7"><span class="toc-item-num">7&nbsp;&nbsp;</span>Julia set</a></span></li><li><span><a href="#Mandelbrot-set" data-toc-modified-id="Mandelbrot-set-8"><span class="toc-item-num">8&nbsp;&nbsp;</span>Mandelbrot set</a></span></li><li><span><a href="#Game-of-Life" data-toc-modified-id="Game-of-Life-9"><span class="toc-item-num">9&nbsp;&nbsp;</span>Game of Life</a></span></li><li><span><a href="#Ising-model-on-the-2-dimensional-square-lattice" data-toc-modified-id="Ising-model-on-the-2-dimensional-square-lattice-10"><span class="toc-item-num">10&nbsp;&nbsp;</span>Ising model on the 2-dimensional square lattice</a></span></li><li><span><a href="#Appendix" data-toc-modified-id="Appendix-11"><span class="toc-item-num">11&nbsp;&nbsp;</span>Appendix</a></span></li></ul></div>

# %%
@time using Distributions
@time using Roots
@time using Plots
@time using Chain
@time using Ising2D
@time using DifferentialEquations
@time using ProgressMeter
VERSION

# %% [markdown]
# ## probability that the confidence interval contains the true value

# %%
using Distributions, Roots, Plots
x ⪅ y = x < y || x ≈ y
pval(d, k) = sum(pdf(d, j) for j in support(d) if pdf(d, j) ⪅ pdf(d, k))
prob_ci_contains_p(d, α=0.05) = sum(pdf(d, k) for k in support(d) if pval(d, k) ≥ α)
n, α, p = 20, 0.05, range(0, 1; length=1001)
plot(p, prob_ci_contains_p.(Binomial.(n, p), α), ylim=(0.94, 1), label="")
plot!(; xlabel="parameter p", xtick=0:0.1:1)
plot!(; ylabel="probability of CI(α=$α) ∋ p")
plot!(; title="true distribution = Binomial($n, p)")
hline!([1 - α]; label="1 - α", legend=:bottomright, ls=:dash)

# %% [markdown]
# ## distribution of confidence intervals of the binomial distribution model

# %%
using Distributions, Roots, Plots
x ⪅ y = x < y || x ≈ y
pval(d, k) = sum(pdf(d, j) for j in support(d) if pdf(d, j) ⪅ pdf(d, k))
ci(n, k, α) = (CI = find_zeros(p -> pval(Binomial(n, p), k) - α, 0, 1);
    isone(length(CI)) ? CI[1] < 0.5 ? (0.0, CI[1]) : (CI[1], 1.0) : (CI[1], CI[end]))
n, p = 20, 0.4
CIs, Freqs, lim = ci.(n, 0:n, 0.05), pdf.(Binomial(n, p), 0:n), (-0.05, 1.05)
plot(; size=(400, 400), xlabel="min of CI", ylabel="max of CI", xlim=lim, ylim=lim)
scatter!(CIs; ms=max.(1, .√(500Freqs)), mz=Freqs, colorbar=false, label="confidence intervals")
plot!([[p] [p]]; c=:red, seriestype=[:hline :vline], label=["true value" ""], legend=:bottomright)

# %% [markdown]
# ## Fisher's exact test version

# %%
using Distributions, Roots, Plots
x ⪅ y = x < y || x ≈ y
pval(d, k) = sum(pdf(d, j) for j in support(d) if pdf(d, j) ⪅ pdf(d, k))
prob_ci_contains_t(d, α=0.05) = sum(pdf(d, k) for k in support(d) if pval(d, k) ≥ α)
ns, nf, n, α, t = 60, 40, 50, 0.05, range(-5, 5; length=1001)
plot(t, @.(prob_ci_contains_t(FisherNoncentralHypergeometric(ns, nf, n, exp(t)), α)), ylim=(0.94, 1), label="")
plot!(; xlabel="parameter t = log ω", xtick=-5:0.5:5)
plot!(; ylabel="probability of CI(α=$α) ∋ t")
plot!(; title="true distribution = FisherNoncentralHypergeometric($ns, $nf, $n, exp(t))", titlefontsize=11)
hline!([1 - α]; label="1 - α", legend=:bottomright, ls=:dash)

# %%
using Distributions, Roots, Plots; x ⪅ y = x < y || x ≈ y
pval(d, k) = sum(pdf(d, j) for j in support(d) if pdf(d, j) ⪅ pdf(d, k))
pval(ns, nf, n, t, k) = pval(FisherNoncentralHypergeometric(ns, nf, n, exp(t)), k)
ci(ns, nf, n, k, α=0.05) = (CI = find_zeros(t -> pval(ns, nf, n, t, k) - α, -20, 20);
    isone(length(CI)) ? CI[1] < 0 ? (-Inf, CI[1]) : (CI[1], Inf) : (CI[1], CI[end]))
d = Hypergeometric(30, 20, 25); k = support(d)[2:end-1]
CIs = ci.(params(d)..., k, 0.05); Freqs = pdf.(d, k)
plot(; size=(400, 400), xlabel="min of CI", ylabel="max of CI", xlim = (-8, 3), ylim = (-3, 8))
scatter!(CIs; ms=max.(1, .√(500Freqs)), mz=Freqs, colorbar=false, label="confidence intervals")
plot!([[0.0] [0.0]]; c=:red, seriestype=[:hline :vline], label=["true value" ""], legend=:topleft)

# %%
using Distributions, Roots, Plots
x ⪅ y = x < y || x ≈ y
pval(d, k) = sum(pdf(d, j) for j in support(d) if pdf(d, j) ⪅ pdf(d, k))
pval(ns, nf, n, t, k) = pval(FisherNoncentralHypergeometric(ns, nf, n, exp(t)), k)
ci(ns, nf, n, k, α=0.05) = (CI = find_zeros(t -> pval(ns, nf, n, t, k) - α, -10, 10); (CI[1], CI[end]))
pval(A::AbstractMatrix) = ((a, b, c, d) = vec(A); pval(a+b, c+d, a+c, 0.0, a))
ci(A::AbstractMatrix, α=0.05) = ((a, b, c, d) = vec(A); exp.(ci(a+b, c+d, a+c, a, α)))
A = [10 10; 7 27]
@show pval(A)
@show ci(A);

# %% [markdown]
# 以下のセルを実行する前に表示されている画像を `R fisher.test.jpg` というファイル名で保存しておくこと.

# %%
using Base64, Markdown
showimg(mime, fn; tag="img") = open(fn) do f
    base64 = base64encode(f)
    display("text/html", """<$tag src="data:$mime;base64,$base64" />""")
end
showimg("image/jpg", "R fisher.test.jpg"; tag="img width=80%")
md"""
R の `fisher.test` ではP値が5%を切っているのに、信頼区間に ``0`` が含まれてしまう場合がある。
上の実装ではそのようなことは起こっていない。
"""

# %% [markdown]
# 上の添付画像の内容はR言語および `RCall.jl` をインストールしておけばJulia内で実行可能である。

# %% [markdown]
# ## 3-dimensional random walk

# %%
using Plots, Chain, Random; Random.seed!(4649373)
N = 3000
X, Y, Z = @chain (cumsum(rand([-1, 1], 3, N), dims=2); @views((_[1,:], _[2,:], _[3,:])))
xlim, ylim, zlim = extrema.((X, Y, Z))
@gif for n in [fill(10, 20); 10:10:N; fill(N, 20)]
    J, K = 1:max(1, n-200), max(1, n-199):n 
    plot(X[K], Y[K], Z[K]; xlim, ylim, zlim, legend=false, color=:red)
    plot!(X[J], Y[J], Z[J]; xlim, ylim, zlim, lw=0.5, alpha=0.5, color=:red)
    title!("3-dimensional random walk")
end

# %% [markdown]
# ## Ising model on the 2-dimensional square lattice

# %%
using Ising2D # Cheat!
gif_mcmc_ising2d()

# %% [markdown]
# ## Kuramoto model

# %%
using Distributions, DifferentialEquations, Plots
kuramoto!(dθ, θ, param, t) = (N = length(θ); (K, ω) = param; 
    for i in 1:N dθ[i] = ω[i] + K*mean(sin(θ[j] - θ[i]) for j in 1:N) end)
(m, n, d, v) = 32, 16, Normal(0, 2), 1.0; (K_c, tmax, a) = 2/π/pdf(d, 0), 25.0, 1.3
θ₀, tspan, K, ω = 2π*rand(m, n), (0.0, tmax), a*K_c, rand(d, m, n) .+ v
sol = solve(ODEProblem(kuramoto!, θ₀, tspan, (K, ω)))
@gif for t in [fill(0, 20); 0:0.1:tmax; fill(tmax, 20)]
    plot(; size=(320, 180), title="Kuramoto model at K = $(a)K_c:   t = $t", titlefontsize=10)
    heatmap!(sin.(sol(t))'; c=:bwr, colorbar=false, frame=false, ticks=false)
end

# %%
using Distributions, DifferentialEquations, Plots
kuramoto!(dθ, θ, param, t) = (N = length(θ); (K, ω) = param; 
    for i in 1:N dθ[i] = ω[i] + K*mean(sin(θ[j] - θ[i]) for j in 1:N) end)
(m, n, d, v) = 32, 16, Normal(0, 2), 1.0; (K_c, tmax, a) = 2/π/pdf(d, 0), 25.0, 0.9
θ₀, tspan, K, ω = 2π*rand(m, n), (0.0, tmax), a*K_c, rand(d, m, n) .+ v
sol = solve(ODEProblem(kuramoto!, θ₀, tspan, (K, ω)))
@gif for t in [fill(0, 20); 0:0.1:tmax; fill(tmax, 20)]
    plot(; size=(320, 180), title="Kuramoto model at K = $(a)K_c:   t = $t", titlefontsize=10)
    heatmap!(sin.(sol(t))'; c=:bwr, colorbar=false, frame=false, ticks=false)
end

# %%
using Distributions, DifferentialEquations, Plots
kuramoto!(dθ, θ, param, t) = (N = length(θ); (K, ω) = param; 
    for i in 1:N dθ[i] = ω[i] + K*mean(sin(θ[j] - θ[i]) for j in 1:N) end)
(m, n, d, v) = 32, 16, Normal(0, 2), 1.0; (K_c, tmax, a) = 2/π/pdf(d, 0), 25.0, 1.0
θ₀, tspan, K, ω = 2π*rand(m, n), (0.0, tmax), a*K_c, rand(d, m, n) .+ v
sol = solve(ODEProblem(kuramoto!, θ₀, tspan, (K, ω)))
@gif for t in [fill(0, 20); 0:0.1:tmax; fill(tmax, 20)]
    plot(; size=(320, 180), title="Kuramoto model at K = $(a)K_c:   t = $t", titlefontsize=10)
    heatmap!(sin.(sol(t))'; c=:bwr, colorbar=false, frame=false, ticks=false)
end

# %%
using Distributions, DifferentialEquations, Plots
kuramoto!(dθ, θ, param, t) = (N = length(θ); (K, ω) = param; 
    for i in 1:N dθ[i] = ω[i] + K*mean(sin(θ[j] - θ[i]) for j in 1:N) end)
(m, n, d, v) = 32, 16, Normal(0, 2), 1.0; (K_c, tmax, a) = 2/π/pdf(d, 0), 25.0, 1.1
θ₀, tspan, K, ω = 2π*rand(m, n), (0.0, tmax), a*K_c, rand(d, m, n) .+ v
sol = solve(ODEProblem(kuramoto!, θ₀, tspan, (K, ω)))
@gif for t in [fill(0, 20); 0:0.1:tmax; fill(tmax, 20)]
    plot(; size=(320, 180), title="Kuramoto model at K = $(a)K_c:   t = $t", titlefontsize=10)
    heatmap!(sin.(sol(t))'; c=:bwr, colorbar=false, frame=false, ticks=false)
end

# %%
using Distributions, DifferentialEquations, Plots
kuramoto!(dθ, θ, param, t) = (N = length(θ); (K, ω) = param; 
    for i in 1:N dθ[i] = ω[i] + K*mean(sin(θ[j] - θ[i]) for j in 1:N) end)
(m, n, d, v) = 32, 16, Normal(0, 2), 1.0; (K_c, tmax, a) = 2/π/pdf(d, 0), 25.0, 1.2
θ₀, tspan, K, ω = 2π*rand(m, n), (0.0, tmax), a*K_c, rand(d, m, n) .+ v
sol = solve(ODEProblem(kuramoto!, θ₀, tspan, (K, ω)))
@gif for t in [fill(0, 20); 0:0.1:tmax; fill(tmax, 20)]
    plot(; size=(320, 180), title="Kuramoto model at K = $(a)K_c:   t = $t", titlefontsize=10)
    heatmap!(sin.(sol(t))'; c=:bwr, colorbar=false, frame=false, ticks=false)
end

# %% [markdown]
# ## Julia set

# %%
using Plots, Chain, ProgressMeter
julia_set_count(c, z, N=255) = (for k in 1:N z = z^2 + c; abs2(z) > 4 && return k end; N)
x, y = range(-1.5, 1.5; length=301), range(-1.05, 1.05; length=231)
ts = @chain (range(0.598, 0.613; length=100); [_; reverse(_[2:end-1])])
prog = Progress(length(ts), 0); colorbar = frame = ticks = false
@gif for t in ts
    jsc = julia_set_count.(complex(-0.38, t), complex.(x', y))
    heatmap(x, y, jsc; c=:gist_earth, size=(480, 336), colorbar, frame, ticks)
    next!(prog)
end

# %% [markdown]
# ## Mandelbrot set

# %%
using Plots, Chain, ProgressMeter
mandelbrot_set_count(c, N=2^10) = (z = zero(c); for k in 1:N z = z^2 + c; abs2(z) > 16 && return k end; N)
ts = @chain (exp.(-7range(0, 1; length=101)); [fill(_[1], 20); _; fill(_[end], 20)])
prog = Progress(length(ts), 0); colorbar = frame = ticks = false
@gif for t in ts
    x, y = @chain (range(-t, t; length=300); (-1.3563 .+ _, 0.0685 .+ _))
    msc = @. log(mandelbrot_set_count(complex(x', y)))
    heatmap(x, y, msc; c=reverse(cgrad(:jet1)), size=(340, 340), colorbar, frame, ticks)
    next!(prog)
end

# %%
using Plots, Chain, ProgressMeter
mandelbrot_set_count(c, N=2^12) = (z = zero(c); for k in 1:N z = z^2 + c; abs2(z) > 16 && return k end; N)
ts = @chain (1.4exp.(-12range(0, 1; length=201)); [fill(_[1], 10); _; fill(_[end], 30)])
prog = Progress(length(ts), 0); colorbar = frame = ticks = false
@gif for t in ts
    x, y = @chain (range(-t, t; length=300); (-0.714684 .+ _, 0.299877 .+ _))
    msc = @. log(mandelbrot_set_count(complex(x', y)))
    heatmap(x, y, msc; c=reverse(cgrad(:jet1)), size=(340, 340), colorbar, frame, ticks)
    next!(prog)
end

# %% [markdown]
# ## Game of Life
#
# See also https://twitter.com/genkuroki/status/1378261871329910786

# %%
using Plots, ProgressMeter, Random; m = n = 200; L = 500; prog = Progress(L, 0); v = zeros(Int8, m, n)
u = zero(v); u[2:end-1, 2:end-1] = rand(Int8[0,1], m-2, n-2)
lifegame!(v, u, m, n) = for j in 2:n-1, i in 2:m-1
    nn = u[i-1,j-1] + u[i,j-1] + u[i+1,j-1] + u[i-1,j] + u[i+1,j] + u[i-1,j+1] + u[i,j+1] + u[i+1,j+1]
    v[i,j] = nn == 3 ? 1 : nn == 2 && !iszero(u[i,j]) ? 1 : 0
end
@gif for t in 1:L
    heatmap(u; size=(440, 440), colorbar=false, frame=false, ticks=false)
    lifegame!(v, u, m, n); u .= v; next!(prog)
end

# %%
using Plots, ProgressMeter, Random; m = n = 200; L = 2000; prog = Progress(L, 0); v = zeros(Int8, m, n)
u = zero(v); u[148:152, 148:152] = [1 1 1 0 1; 1 0 0 0 0; 0 0 0 1 1; 0 1 1 0 1; 1 0 1 0 1]
lifegame!(v, u, m, n) = for j in 2:n-1, i in 2:m-1
    nn = u[i-1,j-1] + u[i,j-1] + u[i+1,j-1] + u[i-1,j] + u[i+1,j] + u[i-1,j+1] + u[i,j+1] + u[i+1,j+1]
    v[i,j] = nn == 3 ? 1 : nn == 2 && !iszero(u[i,j]) ? 1 : 0
end
@gif for t in 1:L
    heatmap(u; size=(440, 440), colorbar=false, frame=false, ticks=false)
    lifegame!(v, u, m, n); u .= v; next!(prog)
end

# %% [markdown]
# ## Ising model on the 2-dimensional square lattice

# %%
using Plots, ProgressMeter, Random; rng = Random.default_rng(); (n, skip, L) = (200, 1, 200) 
Q(i, m) = ifelse(i == 1, m, i-1); P(i, m) = ifelse(i == m, 1, i+1)
ising2d!(s, skip, prob, rng) = ((m, n) = size(s); for _ in 1:skip, j in 1:n, i in 1:m
        sij = s[i,j]; k = sij * (s[Q(i, m), j] + s[P(i, m), j] + s[i, Q(j, n)] + s[i, P(j, n)])
        s[i, j] = ifelse(rand(rng) < prob[k+5], -sij, sij) end)
β = log(1 + √2)/2; prob = ((exp(-2β*k) for k in -4:4)...,); s = rand(Int8[-1,1], n, n); prog = Progress(L, 0)
@gif for t in 0:L
    heatmap(s; size=(240, 260), colorbar=false, frame=false, axis=false, ticks=false)
    title!("$(n)×$(n) × $(skip) skip × $(t)"; titlefontsize=11); ising2d!(s, skip, prob, rng); next!(prog)
end

# %%
using Plots, ProgressMeter, Random; rng = Random.default_rng(); (n, skip, L) = (200, 200, 200) 
Q(i, m) = ifelse(i == 1, m, i-1); P(i, m) = ifelse(i == m, 1, i+1)
ising2d!(s, skip, prob, rng) = ((m, n) = size(s); for _ in 1:skip, j in 1:n, i in 1:m
        sij = s[i,j]; k = sij * (s[Q(i, m), j] + s[P(i, m), j] + s[i, Q(j, n)] + s[i, P(j, n)])
        s[i, j] = ifelse(rand(rng) < prob[k+5], -sij, sij) end)
β = log(1 + √2)/2; prob = ((exp(-2β*k) for k in -4:4)...,); s = rand(Int8[-1,1], n, n); prog = Progress(L, 0)
@gif for t in 0:L
    heatmap(s; size=(240, 260), colorbar=false, frame=false, axis=false, ticks=false)
    title!("$(n)×$(n) × $(skip) skip × $(t)"; titlefontsize=11); ising2d!(s, skip, prob, rng); next!(prog)
end

# %% [markdown]
# ## Appendix
#
# See also
#
# * https://twitter.com/genkuroki/status/1378261871329910786
# * https://twitter.com/genkuroki/status/1383197260821843969
# * https://twitter.com/genkuroki/status/1383217975063285762

# %%
