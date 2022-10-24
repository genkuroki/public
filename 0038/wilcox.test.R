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
set.seed(17)
x = rnorm(20)
y = rnorm(30)
wilcox.test(x, y, conf.int=T)

# %%
median(x)

# %%
median(y)

# %%
median(x) - median(y)

# %%
