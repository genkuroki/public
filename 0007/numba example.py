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
import random

def f():
    if random.random() > 0.5:
        return 1
    else:
        return 1.0

[f() for _ in range(0, 20)]

# %%
import random
import numba

@numba.njit()
def g():
    if random.random() > 0.5:
        return 1
    else:
        return 1.0

[g() for _ in range(0, 20)]

# %%
import random
import numba

@numba.njit()
def h():
    if random.random() > 0.5:
        return 1
    else:
        return 2.0

[h() for _ in range(0, 20)]

# %%
h.inspect_types()

# %%
