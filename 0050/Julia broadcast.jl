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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# * https://x.com/utokyo_bunny/status/1841149355736383937
#   * https://mti-lab.github.io/blog/2024/09/26/broadcast_jp.html

# %% [markdown]
# <img src="IMG_6368.jpeg" width=50%>

# %%
A = [
    1 2
    3 4
    5 6
]

# %%
v = [
    7 8
]

# %%
A .* v

# %% [markdown]
# <img src="IMG_6369.jpeg">

# %%
X = collect(reshape(1:24, 3, 4, 2)) # collectは必要ない

# %%
Y = collect(reshape(1:12, 4, 3)') # collectは必要ない

# %%
X .* Y

# %% [markdown]
# ## 他の例

# %%
f(x, y, z) = x * y + z

# %%
x = [1; 2]

# %%
y = [3;; 4]

# %%
z = [5;;; 6]

# %%
f.(x, y, z)

# %%
