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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# * https://github.com/hyrodium/Desmos.jl
# * https://colab.research.google.com/github/genkuroki/public/blob/main/0055/Desmos.l.ipynb

# %%
haskey(ENV, "COLAB_GPU") && (import Pkg; Pkg.add("Desmos"))
using Desmos

state = @desmos begin
    @text "First example"
    @expression cos(x) color=RGB(0, 0.5, 1) color="#f0f"
    @expression (cosh(t), sinh(t)) parametric_domain=-2..3
end

# %%
