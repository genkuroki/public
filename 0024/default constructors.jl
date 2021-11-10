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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# ## デフォルトで定義されるコンストラクタ

# %%
struct Kuma{T}
    p::Tuple{T, T}
end

# %% [markdown]
# これだけで、
#
# * `Kuma(p::Tuple{T, T}) where T`
# * `Kuma{T}(p) where T`
#
# が定義される。内部コンストラクタを定義するとこれらは定義されなくなる。
#
# 多くの場合に自動的に定義されるデフォルトのコンストラクタで十分である。
#
# 複雑なことをやってくれるコンストラクタは外部コンストラクタで行えばよい。

# %%
methods(Kuma) |> display
methods(Kuma{Int}) |> display

# %% [markdown]
# ## デフォルトの`Kuma{T}(p)`は引数の型変換を自動的にしてくれる。

# %% [markdown]
# `Kuma{T}(p) where T` は `p` の型を自動的に `Tuple{T, T}` 型に変換してくれる。

# %%
x = Kuma{Int}((Float32(2.0), big(3)))

# %%
typeof(x.p)

# %% [markdown]
# ## promotionによる自動型変換に一般化

# %% [markdown]
# `Kuma(p::Tuple{T, U}) where {T, U}` をpromotionによって作成。

# %%
Kuma(p::Tuple{T, U}) where {T, U} = Kuma{promote_type(T, U)}(p)

# %%
methods(Kuma)

# %%
y = Kuma((big(2), 3.0))

# %%
typeof(y.p)

# %% [markdown]
# ## `Kuma(x, y)` を定義

# %% [markdown]
# `Kuma(x, y)` を定義。

# %%
Kuma(x, y) = Kuma((x, y))

# %%
methods(Kuma)

# %%
z = Kuma(big(2), 3.0)

# %%
typeof(z.p)

# %% [markdown]
# ## `==` を定義

# %%
y == z

# %%
Base.:(==)(x::Kuma, y::Kuma) = x.p == y.p

# %%
y == z

# %%
