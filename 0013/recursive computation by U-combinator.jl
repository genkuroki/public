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
Fact = (u -> u(u))(u -> n -> n == 0 ? 1 : n * u(u)(n - 1))
@show [Fact(n) for n in 0:10];

# %%
f(u, n) = n == 0 ? 1 : n * u(u, n - 1)
@show [f(f, n) for n in 0:10];

# %%
@code_warntype f(f, 5)

# %%
@code_typed f(f, 5)

# %%
@code_llvm debuginfo=:none f(f, 5)

# %%
@code_native debuginfo=:none f(f, 5)

# %%
g(u, n, a, b) = n == 0 ? a : u(u, n-1, b, a+b)
@show [g(g, n, 0, 1) for n in 0:10];

# %%
U(u) = u(u)

# %%
U(U)

# %%
