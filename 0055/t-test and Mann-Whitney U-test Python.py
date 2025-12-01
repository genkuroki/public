# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:hydrogen
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# https://okumuralab.org/~okumura/python/t_u_test.html はヒストグラムを使っている点が非常によろしくない。

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

rng = np.random.default_rng()

# %%
M = 10
N = 100000

def p_value():
    data1 = rng.standard_normal(M)
    data2 = rng.standard_normal(M)
    _, p = stats.ttest_ind(data1, data2, equal_var=True)
    return p

%time p_values = [p_value() for _ in range(N)]
plt.hist(p_values, color="lightgray", edgecolor="black",
         bins=np.arange(21)/20, density=True)

# %%
M = 10
N = 100000

def p_value():
    data1 = rng.random(M)
    data2 = rng.random(M)
    _, p = stats.ttest_ind(data1, data2, equal_var=False)
    return p

%time p_values = [p_value() for _ in range(N)]
plt.hist(p_values, color="lightgray", edgecolor="black",
         bins=np.arange(21)/20, density=True)

# %%
M = 10
N = 100000

def p_value():
    data1 = rng.random(M)
    data2 = rng.random(M)
    _, p = stats.mannwhitneyu(data1, data2, method="exact")
    return p

%time p_values = [p_value() for _ in range(N)]
plt.hist(p_values, color="lightgray", edgecolor="black",
         bins=np.arange(21)/20, density=True)

# %%
M = 10
N = 100000

def p_value():
    data1 = rng.random(M)
    data2 = rng.random(M)
    _, p = stats.mannwhitneyu(data1, data2, method="asymptotic")
    return p

%time p_values = [p_value() for _ in range(N)]
plt.hist(p_values, color="lightgray", edgecolor="black",
         bins=np.arange(21)/20, density=True)

# %%
M = 5
N = 100000

def p_value():
    data1 = rng.random(M)
    data2 = rng.random(M)
    _, p = stats.ttest_ind(data1, data2, equal_var=True)
    return p

%time p_values = [p_value() for _ in range(N)]
plt.hist(p_values, color="lightgray", edgecolor="black",
         bins=np.arange(21)/20, density=True)

# %%
M = 3
N = 100000

def p_value():
    data1 = rng.random(M)
    data2 = rng.random(M)
    _, p = stats.ttest_ind(data1, data2, equal_var=True)
    return p

%time p_values = [p_value() for _ in range(N)]
plt.hist(p_values, color="lightgray", edgecolor="black",
         bins=np.arange(21)/20, density=True)

# %%
