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
    display_name: Julia 1.9.2
    language: julia
    name: julia-1.9
---

# 2×2の分割表の独立性に関する一様最強力不偏検定の実装例

* 黒木玄
* 2021-12-24, 2023-08-13

```julia
ENV["COLUMNS"] = 200

using Distributions

using StatsPlots
default(fmt=:png, size=(400, 300), titlefontsize=10, tickfontsize=6)

using SymPy
# Override the Base.show definition of SymPy.jl:
# https://github.com/JuliaPy/SymPy.jl/blob/29c5bfd1d10ac53014fa7fef468bc8deccadc2fc/src/types.jl#L87-L105
@eval SymPy function Base.show(io::IO, ::MIME"text/latex", x::SymbolicObject)
    print(io, as_markdown("\\displaystyle " * sympy.latex(x, mode="plain", fold_short_frac=false)))
end
@eval SymPy function Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{Sym})
    function toeqnarray(x::Vector{Sym})
        a = join(["\\displaystyle " * sympy.latex(x[i]) for i in 1:length(x)], "\\\\")
        """\\left[ \\begin{array}{r}$a\\end{array} \\right]"""
    end
    function toeqnarray(x::AbstractArray{Sym,2})
        sz = size(x)
        a = join([join("\\displaystyle " .* map(sympy.latex, x[i,:]), "&") for i in 1:sz[1]], "\\\\")
        "\\left[ \\begin{array}{" * repeat("r",sz[2]) * "}" * a * "\\end{array}\\right]"
    end
    print(io, as_markdown(toeqnarray(x)))
end

using Memoization

using BenchmarkTools

x ⪅ y = x < y || x ≈ y    # \lessapprox
x ⪆ y = y ⪅ x             # \gtrapprox
x ⪉ y = x < y && !(x ≈ y) # \lnapprox
x ⪊ y = y ⪉ x             # \gnapprox

safemul(x, y) = x == 0 ? zero(float(typeof(x))) : x * y
safediv(x, y) = x == 0 ? zero(float(typeof(x))) : x / y

plot(sin; size=(300, 200))
```

## 2×2の分割表の独立性検定


### 2つの独立な二項分布モデルの場合を考える

```julia
pdfbin(r, s, p, q, x, y) = pdf(Binomial(r, p), x) * pdf(Binomial(s, q), y)
supportbin(r, s) = Iterators.product(0:r, 0:s)
prob_rejection(testfunc, α, r, s, p, q) = 
    sum(pdfbin(r, s, p, q, a, c) * testfunc(a, r-a, c, s-c; α) for (a, c) in supportbin(r, s))
```

### 独立性のχ²検定

```julia
function chisqstat(a, b, c, d)
    s, f, m ,n = a + b, c + d, a + c, b + d
    safediv((a*d - b*c)^2 * (m + n), s * f * m * n)
end
chisqstat(A) = chisqstat(A...)

@syms a b c d
chisqstat(a, b, c, d)
```

```julia
@memoize pval_chisq(a, b, c, d) = ccdf(Chisq(1), chisqstat(a, b, c, d))
pval_chisq(A) = pval_chisq(A...)

test_chisq(a, b, c, d; α = 0.05) = float(pval_chisq(a, b, c, d) < α)
test_chisq(A; α = 0.05) = test_chisq(A...; α)
```

```julia
a, b, c, d = 1, 3, 6, 1
pval_chisq(a, b, c, d), test_chisq(a, b, c, d)
```

### Fisher検定

```julia
@memoize function pval_fisher(a, b, c, d)
    s, f, m ,n = a + b, c + d, a + c, b + d
    (s == 0 || f == 0 || m == 0 || n == 0) && return 1.0
    hg = Hypergeometric(s, f, m)
    p = sum(pdf(hg, j) for j in support(hg) if pdf(hg, j) ⪅ pdf(hg, a))
    min(1, p)
end
pval_fisher(A) = pval_fisher(A...)

test_fisher(a, b, c, d; α = 0.05) = float(pval_fisher(a, b, c, d) < α)
test_fisher(A; α = 0.05) = test_fisher(A...; α)
```

```julia
a, b, c, d = 1, 3, 9, 1
pval_fisher(a, b, c, d), test_fisher(a, b, c, d)
```

```julia
α = 0.05
r, s, p, q = 10, 20, 0.5, 0.5
prob_rejection(test_chisq, α, r, s, p, q), prob_rejection(test_fisher, α, r, s, p, q)
```

### 独立性の一様最強力不偏検定

```julia
function _pvals_umpu(dist, x)
    px = pdf(dist, x)
    P1 = sum(pdf(dist, j) for j in support(dist) if pdf(dist, j) ⪅ px; init = 0.0)
    P0 = sum(pdf(dist, j) for j in support(dist) if pdf(dist, j) ⪉ px; init = 0.0)
    min(1, P0), min(1, P1)
end

function _test_umpu(dist, x; α = 0.05)
    px = pdf(dist, x)
    P1 = sum(pdf(dist, j) for j in support(dist) if pdf(dist, j) ⪅ px; init = 0.0)
    P1 ⪉ α && return 1.0
    P0 = sum(pdf(dist, j) for j in support(dist) if pdf(dist, j) ⪉ px; init = 0.0)
    P0 ⪉ α ? (α - P0) / (P1 - P0) : 0.0
end

@memoize function pvals_umpu(a, b, c, d)
    s, f, m ,n = a + b, c + d, a + c, b + d
    hg = Hypergeometric(s, f, m)
    _pvals_umpu(hg, a)
end
pvals_umpu(A) = pvals_umpu(A...)

@memoize function test_umpu(a, b, c, d; α = 0.05)
    s, f, m ,n = a + b, c + d, a + c, b + d
    hg = Hypergeometric(s, f, m)
    _test_umpu(hg, a; α)
end
test_umpu(A; α = 0.05) = test_umpu(A...; α)
```

```julia
α = 0.05
r, s, p, q = 5, 5, 0.4, 0.4
prob_rejection(test_umpu, α, r, s, p, q)
```

```julia
[(a, r-a, c, s-c) for (a, c) in supportbin(r, s)]
```

```julia
[pdfbin(r, s, p, q, a, c) for (a, c) in supportbin(r, s)]
```

```julia
[test_umpu(a, r-a, c, s-c; α) for (a, c) in supportbin(r, s)]
```

```julia
[pdfbin(r, s, p, q, a, c) * test_umpu(a, r-a, c, s-c; α) for (a, c) in supportbin(r, s)]
```

```julia
α, r, s = 0.1, 5, 10
p = q = 0:0.1:1
prob_rejection.(test_umpu, α, r, s, p', q)
```

```julia
prob_rejection.(test_fisher, α, r, s, p', q)
```

```julia
prob_rejection.(test_chisq, α, r, s, p', q)
```

### 独立性の検定における第一種の過誤の確率のプロット

```julia tags=[]
function plot_realalphas(; r = 10, s = 20, p = 0.4, q = p, kwargs...)
    α = 0:0.005:1
    xtick = ytick = 0:0.1:1
    A = prob_rejection.(test_chisq, α, r, s, p, q)
    B = prob_rejection.(test_fisher, α, r, s, p, q)
    C = prob_rejection.(test_umpu, α, r, s, p, q)    
    P = plot(; legend=:topleft, xtick, ytick)
    plot!(α, A; label="chisq")
    plot!(α, B; label="fisher", ls=:dash)
    plot!(α, C; label="umpu", ls=:dashdot)
    
    α = 0:0.0005:0.1
    xtick = ytick = 0:0.01:1
    A = prob_rejection.(test_chisq, α, r, s, p, q)
    B = prob_rejection.(test_fisher, α, r, s, p, q)
    C = prob_rejection.(test_umpu, α, r, s, p, q)    
    Q = plot(; legend=:topleft, xtick, ytick)
    plot!(α, A; label="chisq")
    plot!(α, B; label="fisher", ls=:dash)
    plot!(α, C; label="umpu", ls=:dashdot)
    
    plot(P, Q; size=(600, 300), layout=(1, 2))
    plot!(; kwargs...)
end
```

```julia tags=[]
plot_realalphas(; r=5, s=5, p=0.3)
```

```julia
plot_realalphas(; r=10, s=10, p=0.3)
```

```julia
plot_realalphas(; r=20, s=20, p=0.3)
```

```julia
plot_realalphas(; r=20, s=50, p=0.3)
```

```julia
plot_realalphas(; r=50, s=50, p=0.3)
```

```julia
plot_realalphas(; r=20, s=100, p=0.3)
```

```julia
plot_realalphas(; r=50, s=100, p=0.3)
```

### 独立性の検定における検出力のプロット

```julia
function plot_powers(; α = 0.05, r = 10, s = 20,
        levels = [0.0:0.4α:2α; 0.2:0.2:1], kwargs...)
    p = q = 0:0.01:1
    xtick = ytick = 0:0.1:1
    
    C = prob_rejection.(test_umpu, α, r, s, p', q)
    A = prob_rejection.(test_chisq, α, r, s, p', q)
    B = prob_rejection.(test_fisher, α, r, s, p', q)
    
    R = contour(p, q, C; clim=(0, 1), title="power of umpu", levels, xtick, ytick)
    P = contour(p, q, A; clim=(0, 1), title="power of chisq", levels, xtick, ytick)
    Q = contour(p, q, B; clim=(0, 1), title="power of fisher", levels, xtick, ytick)
    
    AC = @. log2(C) - log2(A)
    AB = @. log2(B) - log2(A)
    CB = @. log2(B) - log2(C)
    
    PR = heatmap(p, q, AC; clim=(-1, 1), c=:bwr, title="power: chisq b-w-r umpu", xtick, ytick)
    PQ = heatmap(p, q, AB; clim=(-1, 1), c=:bwr, title="power: chisq b-w-r fisher", xtick, ytick)
    RQ = heatmap(p, q, CB; clim=(-1, 1), c=:bwr, title="power: umpu b-w-r fisher", xtick, ytick)

    plot(R, P, Q, PR, PQ, RQ; size=(800, 560), layout=(2, 3),
        colorbar=false, contour_labels=true)
    plot!(; kwargs...)
end
```

```julia
plot_powers(; r = 5, s = 5)
```

```julia
plot_powers(; r = 10, s = 10)
```

```julia
plot_powers(; r = 20, s = 20)
```

```julia
plot_powers(; r = 20, s = 50)
```

```julia
plot_powers(; r = 50, s = 50)
```

```julia
plot_powers(; r = 20, s = 100)
```

```julia
plot_powers(; r = 50, s = 100)
```

```julia

```
