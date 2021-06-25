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
#     display_name: Julia 1.7.0-beta2
#     language: julia
#     name: julia-1.7
# ---

# %%
VERSION

# %%
f((a, b, c)) = @show a b c
f(1:10);

# %%
g((; b, d, a)) = @show a b d
g((a = 1, b = 2, c = 3, d = 4, e = 5));

# %%
h((a, ((b, c), d))) = @show a b c d
h((1, ((2, 3, 4), 5, 6), 7, 8));

# %%
k(((; b, a), ((c, d), (; f, e)))) = @show a b c d e f
k(((a=1, b=2, g=7), ((3, 4, 8), (e=5, f=6, h=9), 10)));

# %%
(((x, y),) -> y - x^2)((5, 12))

# %%
:((x, y) -> y - x^2) |> Base.remove_linenums!

# %%
:(((x, y)) -> y - x^2) |> Base.remove_linenums!

# %%
:(((x, y),) -> y^2 - x) |> Base.remove_linenums!

# %%
using OffsetArrays

# %%
v = OffsetArray(-3:3, -3:3)

# %%
a, b, c = v
a, b, c

# %%
(((x, y),) -> y - x^2)(v)

# %%
