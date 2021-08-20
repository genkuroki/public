# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
]activate .

# %%
]st

# %%
]add CairoMakie LaTeXStrings

# %%
]st

# %%
using CairoMakie, LaTeXStrings

xs = 0.0:0.01:2.0
lines(xs, sin.(π.*xs), axis=(title=L"\sin{x}",))

# %%
using Pkg
println("Julia v", VERSION)
Pkg.status("CairoMakie")
Pkg.status("Makie"; mode=PKGMODE_MANIFEST)
Pkg.status("LaTeXStrings")

# %%
