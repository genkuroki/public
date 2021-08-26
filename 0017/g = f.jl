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
using BenchmarkTools

f(x) = x^2
g = f
const g_const = f

F(x) = f(x)
G(x) = g(x)
G_const(x) = g_const(x)

@btime f(2.1)
@btime F(2.1)
@btime G(2.1)
@btime G_const(2.1)
println()
@btime g(2.1)
@btime $g(2.1)
println()
@btime F($(Ref(2.1))[])
@btime G($(Ref(2.1))[])
@btime G_const($(Ref(2.1))[]);

# %%
@code_llvm debuginfo=:none f(2.1)

# %%
@code_llvm debuginfo=:none F(2.1)

# %%
@code_llvm debuginfo=:none G(2.1)

# %%
@code_llvm debuginfo=:none G_const(2.1)

# %%
