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

# %%
import numpy as np

# %%
# A.shape == (3, 2)
A = np.array([[1, 2],
              [3, 4],
              [5, 6]])
A

# %%
# v.shape == (1, 2)
v = np.array([[7, 8]])
v

# %%
# Broadcast & multiply
A * v
# > array([[ 7, 16],
# >        [21, 32],
# >        [35, 48]])

# %%
X = np.array([[[1, 4, 7, 10],
               [2, 5, 8, 11],
               [3, 6, 9, 12]],
              [[13, 16, 19, 22],
               [14, 17, 20, 23],
               [15, 18, 21, 24]]])
X

# %%
Y = np.array([[1, 2, 3, 4],
           [5, 6, 7, 8],
           [9, 10, 11, 12]])
Y

# %%
X * Y

# %%
