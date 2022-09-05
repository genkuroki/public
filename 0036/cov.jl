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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown] tags=[]
# 次の公式の確認: $\bar{x} = (x_1+\cdots+x_n)/n$, $\bar{y} = (y_1+\cdots+y_n)/n$ とおくと,
#
# $$
# \frac{1}{n(n-1)}\sum_{1\le i<j\le n} (x_i - x_j)(y_i - y_j) =
# \frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y}).
# $$

# %%
using Distributions

# %%
function mycov(x, y)
    @assert length(x) == length(y)
    n = length(x)
    sum((x[i] - x[j])*(y[i] - y[j]) for i in 1:n for j in i+1:n) / (n*(n-1))
end

# %%
x = [1, 2, 3, 4]
y = [4, 2, -1, 3];

# %%
cov(x, y)

# %%
mycov(x, y)

# %%
X, Y = randn(10), randn(10);

# %%
cov(X, Y)

# %%
mycov(X, Y)

# %%
all((X = randn(10); Y = randn(10); cov(X, Y) ≈ mycov(X, Y)) for _ in 1:10^4)

# %%
using SymPy
n = 4
x, y = @eval @syms x[1:$n]::real y[1:$n]::real

# %%
cov_xy = cov(x, y).simplify()

# %%
mycov_xy = mycov(x, y).expand()

# %%
cov_xy == mycov_xy

# %%
n = 10
x, y = @eval @syms x[1:$n]::real y[1:$n]::real
cov(x, y).simplify() == mycov(x, y).expand()

# %%
