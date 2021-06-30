# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:hydrogen
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
tak(x, y, z) = x ≤ y ? z : tak(tak(x - 1, y, z), tak(y - 1, z, x), tak(z - 1, x, y))

print("Julia ", VERSION, ':')
@time tak(20, 10, 0)
print("Julia ", VERSION, ':')
@time tak(20, 10, 0)
print("Julia ", VERSION, ':')
@time tak(20, 10, 0)

# %%
using BenchmarkTools
print("Julia ", VERSION, ':')
@btime tak(20, 10, 0)

# %%
from time import time

def tak(x, y, z):
    if x <= y: return z
    return tak(tak(x - 1, y, z), tak(y - 1, z, x), tak(z - 1, x, y))

start = time(); res = tak(20, 10, 0); time() - start

# %%
start = time(); res = tak(20, 10, 0); time() - start

# %%
start = time(); res = tak(20, 10, 0); time() - start

# %%
from time import time
from numba import jit

@jit(nopython=True)
def tak(x, y, z):
    if x <= y: return z
    return tak(tak(x - 1, y, z), tak(y - 1, z, x), tak(z - 1, x, y))

start = time(); res = tak(20, 10, 0); time() - start

# %%
start = time(); res = tak(20, 10, 0); time() - start

# %%
start = time(); res = tak(20, 10, 0); time() - start

# %%
