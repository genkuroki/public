# -*- coding: utf-8 -*-
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
.Table <- matrix(c(15, 13, 5, 15), 2, 2, byrow=TRUE)
print(.Table)

# epiR で分析
analysis <- epiR::epi.2by2(
  dat = .Table,
  method = "cross.sectional",
  digits = 3,
  conf.level = 0.95,
  units = 1,
  outcome = "as.columns"
)

# Wald法
rr_result <- analysis$massoc.detail$PR.strata.wald
cat("\nWald RR:\n")
print(rr_result)

# スコア法
rr_result <- analysis$massoc.detail$PR.strata.score
cat("\nScore RR:\n")
print(rr_result)

# %%
(15/(15+13))/(5/(5+15))

# %%
(15/(15+5))/(13/(13+15))

# %%
.Table <- matrix(c(15, 13, 5, 15), 2, 2, byrow=TRUE)

# epiR で分析
analysis <- epiR::epi.2by2(
  dat = .Table,
  method = "cross.sectional",
  digits = 3,
  conf.level = 0.95,
  units = 1,
  outcome = "as.columns"
)
print(analysis)

# Wald法
rr_result <- analysis$massoc.detail$PR.strata.wald
cat("\nWald RR:\n")
print(rr_result)

# スコア法
rr_result <- analysis$massoc.detail$PR.strata.score
cat("\nScore RR:\n")
print(rr_result)

# %%
str(analysis)

# %%
.Table <- matrix(c(15, 13, 5, 15), 2, 2, byrow=FALSE)

# epiR で分析
analysis <- epiR::epi.2by2(
  dat = .Table,
  method = "cross.sectional",
  digits = 3,
  conf.level = 0.95,
  units = 1,
  outcome = "as.columns"
)
print(analysis)

# Wald法
rr_result <- analysis$massoc.detail$PR.strata.wald
cat("\nWald RR:\n")
print(rr_result)

# スコア法
rr_result <- analysis$massoc.detail$PR.strata.score
cat("\nScore RR:\n")
print(rr_result)

# %%
str(analysis)

# %%
