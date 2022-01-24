---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.3
  kernelspec:
    display_name: Julia 1.7.1
    language: julia
    name: julia-1.7
---

# 一様最強力不偏検定の実装

* 黒木玄
* 2021-01-24

$P(x|\theta)$ は二項分布またはFisherの非心超幾何分布などの指数型離散分布族の確率函数であるとする:

$$
P(x|\theta) = c(\theta) h(x) e^{f(\theta)x}.
$$

ここで, $f(\theta)$ は $\theta$ の狭義単調増加函数であるとする. 例えば, 二項分布の場合には

$$
P(x|\theta) = \binom{n}{x} \theta^x (1 - \theta)^{n-x}
$$

なので

$$
f(\theta) = \log\frac{\theta}{1 - \theta}, \quad
c(\theta) = (1 - \theta)^n, \quad
h(x) = \binom{n}{x}
$$


このとき, $1 = \sum_x P(x|\theta)$ の両辺を $\theta$ で偏微分すると,

$$
\frac{P_\theta(x|\theta)}{P(x|\theta)} =
\frac{\partial}{\partial\theta}\log P(x|\theta) =
\frac{c'(\theta)}{c(\theta)} + f'(\theta)x
$$

より, $X \sim P(\cdot|\theta)$ のとき,

$$
0 =
\sum_x \left(\frac{c'(\theta)}{c(\theta)} + f'(\theta)x\right)P(x|\theta) =
\frac{c'(\theta)}{c(\theta)} + f'(\theta)E[X|\theta].
$$

ゆえに,

$$
P_\theta(x|\theta) = f'(\theta)(x - E[X|\theta]) P(x|\theta).
$$


帰無仮説 $\theta = \theta_0$, 対立仮説 $\theta \ne \theta_0$ の一様最強力不偏検定の検定函数(データ $x$ に対して帰無仮説 $\theta = \theta_0$ の棄却確率を返す函数)は以下の形で与えられることが知られている:

$$
\phi(x|\theta_0) =
\begin{cases}
1 & (x < a\; \text{or}\; b < x) \\
r_a & (x = a) \\
r_b & (x = b) \\
0 & (a < x < b) \\
\end{cases}
$$

ここで $a < b$ と $0\le r_a, r_b\le 1$ は次の条件で特徴付けられる:

$$
E[\phi(X|\theta_0)|\theta_0] = \alpha, \quad
\left.\frac{\partial}{\partial\theta}\right|_{\theta=\theta_0} E[\phi(X|\theta_0)|\theta] = 0.
$$

これらの等式は以下のように書き下される:

$$
\begin{aligned}
&
P(a|\theta_0)r_a + P(b|\theta_0)r_b + \sum_{x < a\; \text{or}\; b < x} P(x|\theta_0) = \alpha,
\\ &
\left.\frac{\partial}{\partial\theta}\right|_{\theta=\theta_0}
\left(P(a|\theta)r_a + P(b|\theta)r_b + \sum_{x < a\; \text{or}\; b < x} P(x|\theta)\right) = 0.
\end{aligned}
$$


後者の等式は, 上で示したように　$P_\theta(x|\theta) = f'(\theta)(x - E[X|\theta]) P(x|\theta)$ なので, 次と同値になることがわかる:

$$
a P(a|\theta_0) r_a + b P(b|\theta_0) r_b + \sum_{x < a\; \text{or}\; b < x} x P(x|\theta_0) =
E[X|\theta_0]E[\phi(X|\theta_0)|\theta_0] = \alpha E[X|\theta_0].
$$


したがって, $r_a, r_b$ は次の連立一次方程式の解になる:

$$
\begin{aligned}
&
P(a|\theta_0)r_a + P(b|\theta_0)r_b =
\alpha - \sum_{x < a\; \text{or}\; b < x} P(x|\theta_0) =: y_1,
\\ &
a P(a|\theta_0) r_a + b P(b|\theta_0) r_b =
\alpha E[X|\theta_0] - \sum_{x < a\; \text{or}\; b < x} x P(x|\theta_0) =: y_2.
\end{aligned}
$$


すなわち, 

$$
r_a = \frac{P(b|\theta_0)(b y_1 - y_2)}{(b - a)P(a|\theta_0)P(b|\theta_0)}, \quad
r_b = \frac{P(a|\theta_0)(y_2 - a y_1)}{(b - a)P(a|\theta_0)P(b|\theta_0)}.
$$

これらが $0$ 以上 $1$ 以下になる $a < b$ を探せばよい.

以下のコードはこれをそのまま素直に実装したものになっている.

```julia
ENV["LINES"] = 200
ENV["COLUMNS"] = 200

using Distributions
using Memoization
using RCall
using StatsPlots
plot(Normal(); size=(300, 200))
```

```julia
x ⪅ y = x < y || x ≈ y
x ⪆ y = x > y || x ≈ y

@memoize function pval_exact(dist, k)
    pk = pdf(dist, k)
    P = sum(pdf(dist, j) for j in support(dist) if pdf(dist, j) ⪅ pk; init = 0.0)
    min(1, P)
end

@memoize pval_eqtailed(dist, k) = min(1, 2cdf(dist, k), 2ccdf(dist, k-1))

@memoize function pval_normal(dist, k)
    μ, σ = mean(dist), std(dist)
    2ccdf(Normal(), abs(k - μ)/σ)
end
```

```julia
function solve_umpu(dist, α = 0.05)
    s = support(dist)
    pmf(x) = pdf(dist, x)
    sol = Vector{Float64}[]
    r = Vector{Float64}(undef, length(s))
    for i in 1:length(s), j in length(s):-1:i+1
        a, b = s[i], s[j]
        P0 = cdf(dist, a - 1) + ccdf(dist, b)
        P0 > α && continue
        M0 = sum(s[k]*pmf(s[k]) for k in 1:i-1; init=0.0) + sum(s[k]*pmf(s[k]) for k in j+1:length(s); init=0.0)
        y1 = α - P0
        y2 = mean(dist)*α - M0
        d = (b - a)*pmf(a)*pmf(b)
        ri = pmf(b)*(b*y1 - y2)/d
        rj = pmf(a)*(y2 - a*y1)/d
        !(0 ⪅ ri ⪅ 1 && 0 ⪅ rj ⪅ 1) && continue 
        r .= 0
        r[1:i-1] .= 1
        r[j+1:end] .= 1
        r[i] = ri
        r[j] = rj
        push!(sol, copy(r))
    end
    sol
end
```

```julia
n = 10
theta = 0.1
dist = Binomial(n, theta)
alpha = 0.05
sol = solve_umpu(dist, alpha)
```

```julia
R"library(ump)"
@rput n theta alpha
r_R = Float64[]
for k in 0:n
    @rput k
    R"rk = umpu.binom(k, n, theta, alpha)"
    @rget rk
    push!(r_R, rk)
end
[r_R]
```

```julia
[r_R - sol[1]]
```

```julia
@memoize function umpu!(r, dist, α = 0.05)
    s = support(dist)
    pmf(x) = pdf(dist, x)
    r .= 0
    for i in eachindex(s), j in lastindex(s):-1:i+1
        a, b = s[i], s[j]
        P0 = cdf(dist, a - 1) + ccdf(dist, b)
        P0 > α && continue
        M0 = sum(s[k]*pmf(s[k]) for k in 1:i-1; init=0.0) + sum(s[k]*pmf(s[k]) for k in j+1:length(s); init=0.0)
        y1 = α - P0
        y2 = mean(dist)*α - M0
        d = (b - a)*pmf(a)*pmf(b)
        ri = pmf(b)*(b*y1 - y2)/d
        rj = pmf(a)*(y2 - a*y1)/d
        !(0 ⪅ ri ⪅ 1 && 0 ⪅ rj ⪅ 1) && continue 
        r[1:i-1] .= 1
        r[j+1:end] .= 1
        r[i] = ri
        r[j] = rj
        break
    end
    r
end

@memoize function umpu(dist, α = 0.05)
    r = umpu!(zeros(length(support(dist))), dist, α)
    s = support(dist)
    function phi(x)
        i = x - minimum(s) + 1
        1 ≤ i ≤ lastindex(s) && return r[i]
        1.0
    end
    phi
end

function verify_umpu(dist, α = 0.05)
    phi = umpu(dist, α)
    s = support(dist)
    phix = phi.(s)
    EphiX = sum(x -> phi(x)*pdf(dist, x), s)
    EXphiX = sum(x -> x*phi(x)*pdf(dist, x), s)
    
    println("0 ≤ ϕ(x) ≤ 1     is ", all(0 .≤ phix .≤ 1))
    println("E[ϕ(X)]  ≈ α     is ", EphiX ≈ α)
    println("E[Xϕ(X)] ≈ αE[X] is ", EXphiX ≈ α*mean(dist))
    println()
    println("ϕ(x)          = ", phix)
    println("ϕ_normal(x)   = ", float(pval_normal.(dist, support(dist)) .< α))
    println("ϕ_exact(x)    = ", float(pval_exact.(dist, support(dist)) .< α))
    println("ϕ_eqtailed(x) = ", float(pval_eqtailed.(dist, support(dist)) .< α))
    println()
    println("ϕ(x) - ϕ_normal(x)   = ", phix - (pval_normal.(dist, s) .< α))
    println("ϕ(x) - ϕ_exact(x)    = ", phix - (pval_exact.(dist, s) .< α))
    println("ϕ(x) - ϕ_eqtailed(x) = ", phix - (pval_eqtailed.(dist, s) .< α))
end
```

竹内彰通『現代数理統計学』問8.9

* https://twitter.com/arts_lib/status/1475435871134707714

```julia
dist = Binomial(4, 1/3)
α = 0.1
verify_umpu(dist, α)
```

```julia
dist = Binomial(4, 0.2)
α = 0.1
verify_umpu(dist, α)
```

```julia
dist = Binomial(4, 0.2)
α = 0.1
pdf.(dist, 0:4)
```

```julia
@show a = 0
@show b = 2
@show P_left = sum(pdf.(dist, 0:a-1))
@show P_right = sum(pdf.(dist, b+1:n))
@show P0 = P_left + P_right
@show M0 = sum(x * pdf(dist, x) for x in 3:4)
@show y1 = α - P0
@show y2 = mean(dist)*α - M0
@show d = (b - a)*pdf(dist, a)*pdf(dist, b)
@show r_a = pdf(dist, b)*(b*y1 - y2)/d
@show r_b = pdf(dist, a)*(y2 - a*y1)/d;
```

```julia
@show a = 0
@show b = 3
@show P_left = sum(pdf.(dist, 0:a-1))
@show P_right = sum(pdf.(dist, b+1:n))
@show P0 = P_left + P_right
@show M0 = sum(x * pdf(dist, x) for x in 0:a-1; init=0.0) + sum(x * pdf(dist, x) for x in b+1:n; init=0.0)
@show y1 = α - P0
@show y2 = mean(dist)*α - M0
@show d = (b - a)*pdf(dist, a)*pdf(dist, b)
@show r_a = pdf(dist, b)*(b*y1 - y2)/d
@show r_b = pdf(dist, a)*(y2 - a*y1)/d;
```

```julia
dist = Binomial(10, 0.05)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.1)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.15)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.2)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.25)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.3)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.35)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.4)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.45)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.5)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = Binomial(10, 0.55)
α = 0.05
verify_umpu(dist, α)
```

```julia
dist = FisherNoncentralHypergeometric(11, 12, 15, 0.9)
α = 0.05
verify_umpu(dist, α)
```

```julia
prob_rejection(testfunc, α, n, p, q = p) = sum(testfunc(n, k, p, α) * pdf(Binomial(n, q), k) for k in 0:n)

function plot_bintest_powers(;
        α = 0.05,
        n = 10,
        test1 = (n, k, p, α) -> pval_exact(Binomial(n, p), k) < α,    title1 = "exact",
        test2 = (n, k, p, α) -> pval_eqtailed(Binomial(n, p), k) < α, title2 = "eqtailed",
        levels = [0.0:0.4α:2α; 0.2:0.2:1], kwargs...)
    p = q = 0.01:0.01:0.99
    xtick = ytick = 0:0.1:1
    xlim = ylim = (0, 1)
    
    A = prob_rejection.(test1, α, n, p', q)
    B = prob_rejection.(test2, α, n, p', q)
    AB = @. log2(B) - log2(A)
    
    P = contourf(p, q, A; clim=(0, 1), title=title1, levels, xtick, ytick, xlim, ylim, c=:CMRmap)
    Q = contourf(p, q, B; clim=(0, 1), title=title2, levels, xtick, ytick, xlim, ylim, c=:CMRmap)
    PQ = heatmap(p, q, AB; clim=(-1, 1), c=:bwr, title="power: $title1 b-w-r $title2", xtick, ytick, xlim, ylim)
    T1Err = plot(;xtick, ytick=0:0.2α:2α, xlim, ylim=(0, 2α), legend=:top)
    plot!(p, [A[i, i] for i in eachindex(p)]; label=title1)
    plot!(p, [B[i, i] for i in eachindex(p)]; label=title2)
    title!("probability of type-I errors")
    
    plot(P, Q, PQ, T1Err; size=(640, 640), layout=(2, 2), colorbar=false, contour_labels=true)
    plot!(; tickfontsize=5, titlefontsize=8, legendfontsize=7)
end
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 10,
    test1 = (n, k, p, α) -> pval_normal(Binomial(n, p), k) < α, title1 = "normal",
    test2 = (n, k, p, α) -> pval_exact(Binomial(n, p), k) < α,  title2 = "exact",
)
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 10,
    test1 = (n, k, p, α) -> pval_normal(Binomial(n, p), k) < α,   title1 = "normal",
    test2 = (n, k, p, α) -> pval_eqtailed(Binomial(n, p), k) < α, title2 = "eqtailed",
)
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 10,
    test1 = (n, k, p, α) -> pval_exact(Binomial(n, p), k) < α,    title1 = "exact", 
    test2 = (n, k, p, α) -> pval_eqtailed(Binomial(n, p), k) < α, title2 = "eqtailed",
)
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 10,
    test1 = (n, k, p, α) -> umpu(Binomial(n, p), α)(k),         title1 = "umpu",
    test2 = (n, k, p, α) -> pval_normal(Binomial(n, p), k) < α, title2 = "normal",
)
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 10,
    test1 = (n, k, p, α) -> umpu(Binomial(n, p), α)(k),         title1 = "umpu",
    test2 = (n, k, p, α) -> pval_exact(Binomial(n, p), k) < α,  title2 = "exact",
)
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 10,
    test1 = (n, k, p, α) -> umpu(Binomial(n, p), α)(k),           title1 = "umpu",
    test2 = (n, k, p, α) -> pval_eqtailed(Binomial(n, p), k) < α, title2 = "eqtailed",
)
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 100,
    test1 = (n, k, p, α) -> umpu(Binomial(n, p), α)(k),         title1 = "umpu",
    test2 = (n, k, p, α) -> pval_normal(Binomial(n, p), k) < α, title2 = "normal",
)
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 100,
    test1 = (n, k, p, α) -> umpu(Binomial(n, p), α)(k),         title1 = "umpu",
    test2 = (n, k, p, α) -> pval_exact(Binomial(n, p), k) < α,  title2 = "exact",
)
```

```julia
plot_bintest_powers(
    α = 0.05,
    n = 100,
    test1 = (n, k, p, α) -> umpu(Binomial(n, p), α)(k),           title1 = "umpu",
    test2 = (n, k, p, α) -> pval_eqtailed(Binomial(n, p), k) < α, title2 = "eqtailed",
)
```

```julia
function plot_bintest_phi(ϕ, n; kwargs...)
    x = 0:n
    θ = 0.001:0.001:0.999
    z = ϕ.(x', n, θ)
    heatmap(x, θ, z; ylim=(0, 1), kwargs...)
end
```

```julia
n = 10
α = 0.05

P1 = plot_bintest_phi((x, n, θ) -> umpu(Binomial(n, θ), α)(x), n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="umpu")

P2 = plot_bintest_phi((x, n, θ) -> pval_normal(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="normal")

P3 = plot_bintest_phi((x, n, θ) -> pval_exact(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="exact")

P4 = plot_bintest_phi((x, n, θ) -> pval_eqtailed(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="eqtailed")

plot(P1, P2, P3, P4; size=(640, 640))
plot!(; tickfontsize=7, titlefontsize=10)
```

```julia
n = 10
α = 0.05

P1 = plot_bintest_phi((x, n, θ) -> 1 - umpu(Binomial(n, θ), α)(x), n; clim=(0, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="umpu")

P2 = plot_bintest_phi((x, n, θ) -> pval_normal(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="normal")

P3 = plot_bintest_phi((x, n, θ) -> pval_exact(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="exact")

P4 = plot_bintest_phi((x, n, θ) -> pval_eqtailed(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="eqtailed")

plot(P1, P2, P3, P4; size=(640, 640))
plot!(; tickfontsize=7, titlefontsize=10)
```

```julia
n = 20
α = 0.05

P1 = plot_bintest_phi((x, n, θ) -> umpu(Binomial(n, θ), α)(x), n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="umpu")

P2 = plot_bintest_phi((x, n, θ) -> pval_normal(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="normal")

P3 = plot_bintest_phi((x, n, θ) -> pval_exact(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="exact")

P4 = plot_bintest_phi((x, n, θ) -> pval_eqtailed(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="eqtailed")

plot(P1, P2, P3, P4; size=(640, 640))
plot!(; tickfontsize=5, titlefontsize=10)
```

```julia
n = 20
α = 0.05

P1 = plot_bintest_phi((x, n, θ) -> 1 - umpu(Binomial(n, θ), α)(x), n; clim=(0, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="umpu")

P2 = plot_bintest_phi((x, n, θ) -> pval_normal(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="normal")

P3 = plot_bintest_phi((x, n, θ) -> pval_exact(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="exact")

P4 = plot_bintest_phi((x, n, θ) -> pval_eqtailed(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="eqtailed")

plot(P1, P2, P3, P4; size=(640, 640))
plot!(; tickfontsize=5, titlefontsize=10)
```

```julia
n = 100
α = 0.05

P1 = plot_bintest_phi((x, n, θ) -> umpu(Binomial(n, θ), α)(x), n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:10:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="umpu")

P2 = plot_bintest_phi((x, n, θ) -> pval_normal(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:10:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="normal")

P3 = plot_bintest_phi((x, n, θ) -> pval_exact(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:10:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="exact")

P4 = plot_bintest_phi((x, n, θ) -> pval_eqtailed(Binomial(n, θ), x) < α, n; 
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:10:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="eqtailed")

plot(P1, P2, P3, P4; size=(640, 640))
plot!(; tickfontsize=5, titlefontsize=10)
```

```julia
n = 100
α = 0.05

P1 = plot_bintest_phi((x, n, θ) -> 1 - umpu(Binomial(n, θ), α)(x), n; clim=(0, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:10:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="umpu")

P2 = plot_bintest_phi((x, n, θ) -> pval_normal(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:10:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="normal")

P3 = plot_bintest_phi((x, n, θ) -> pval_exact(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:10:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="exact")

P4 = plot_bintest_phi((x, n, θ) -> pval_eqtailed(Binomial(n, θ), x), n; clim=(-0.1, 1),
    size=(400, 400), colorbar=false, c=:CMRmap, xtick=0:10:n, ytick=0:0.1:1,
    xlabel="x", ylabel="θ", title="eqtailed")

plot(P1, P2, P3, P4; size=(640, 640))
plot!(; tickfontsize=5, titlefontsize=10)
```

```julia

```
