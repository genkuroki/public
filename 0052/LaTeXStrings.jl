# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.11.3
#     language: julia
#     name: julia-1.11
# ---

# %%
using LaTeXStrings

# %%
L"foo $\alpha$" |> display
L"foo $\alpha$" |> dump

# %%
L"foo" |> display
L"foo" |> dump

# %%
LaTeXString("foo") |> display
LaTeXString("foo") |> dump

# %%
L"\pi = %$(float(π))" |> display
L"\pi = %$(float(π))" |> dump

# %%
using Plots
pgfplotsx()
Plots.pgfx_sanitize_string(L"5.5%") |> display
Plots.pgfx_sanitize_string(L"5.5%") |> dump
Plots.pgfx_sanitize_string(L"5.5%") == L"5.5\%"

# %%
