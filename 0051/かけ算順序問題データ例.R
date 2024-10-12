# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:hydrogen
#     text_representation:
#       extension: .R
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown]
# https://x.com/taifu21/status/1588097591606415360
#
# <img src="IMG_6573.png">

# %%
library(exact2x2)
A = matrix(c(10, 7, 21, 47), nrow=2)
A

# %%
exact2x2(A, tsmethod="central", plot=T)

# %%
exact2x2(A, tsmethod="minlik", plot=T)

# %%
