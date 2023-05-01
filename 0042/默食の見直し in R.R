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
# * https://twitter.com/eggplantmed/status/1629147688053923840
#   * https://www.pref.chiba.lg.jp/kyouiku/anzen/hokenn/covid-19.html
#     * https://www.pref.chiba.lg.jp/kyouiku/anzen/hokenn/documents/mokushoku-minaoshi.pdf

# %%
library(exact2x2)

# %%
A = matrix(c(5, 551-5, 27, 1950-27), byrow=T, nrow=2)
A

# %%
chisq.test(A, correct=F)

# %%
exact2x2(A, plot=T)

# %%
5*A

# %%
chisq.test(5*A, correct=F)

# %%
exact2x2(5*A, plot=T)

# %%
6*A

# %%
chisq.test(6*A, correct=F)

# %%
exact2x2(6*A, plot=T)

# %%
10*A

# %%
chisq.test(10*A, correct=F)

# %%
exact2x2(10*A, plot=T)

# %%
B = matrix(c(7*5, 7*(551-5), 2*27, 2*(1950-27)), byrow=T, nrow=2)
B

# %%
chisq.test(B, correct=F)

# %%
exact2x2(B, plot=T)

# %%
2*B

# %%
chisq.test(2*B, correct=F)

# %%
exact2x2(2*B, plot=T)

# %%
