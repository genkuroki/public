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

# %%
N = []
for x = 1:2026, y = x:2026, z = y:2026
    if 2026^2 == x^2+y^2+z^2 push!(N,[x,y,z]) end
end
N

# %%
NN = [[x,y,z] for x in 1:2026 for y in x:2026 for z in y:2026 if 2026^2 == x^2 + y^2 + z^2]

# %%
NN == N

# %%
