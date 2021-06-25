# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# https://docs.julialang.org/en/v1/manual/running-external-programs/

# %%
pr(x) = print(read(x, String))

# %%
pr(`ls -l`)

# %%
ENV["HOGE"] = "hoge"

# %%
pr(pipeline(`env`, `grep '^H'`))

# %%
