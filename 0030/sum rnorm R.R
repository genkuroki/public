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

# %%
version

# %%
RNGkind()

# %%
RNGkind(normal.kind = "Inversion") # default
system.time(sum(rnorm(10^8)))

# %%
RNGkind(normal.kind = "Box-Muller")
system.time(sum(rnorm(10^8)))

# %%
RNGkind(normal.kind = "Ahrens-Dieter")
system.time(sum(rnorm(10^8)))

# %%
RNGkind(normal.kind = "Buggy Kinderman-Ramage")
system.time(sum(rnorm(10^8)))

# %%
RNGkind(normal.kind = "Kinderman-Ramage")
system.time(sum(rnorm(10^8)))

# %%
?RNGkind

# %%
library(RcppZiggurat)

# %%
system.time(sum(zrnorm(10^8)))

# %%
gc()

# %%
system.time(sum(zrnorm(10^8)))

# %%
gc()

# %%
system.time(sum(zrnorm(10^8)))

# %%
?zrnorm

# %%
