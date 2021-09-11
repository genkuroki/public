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
using QuadGK

# %%
a = 1.0
f = x -> sin(a*x)/(a*x)
@time quadgk(f, -1e8, 0, 1e8)
@time quadgk(f, -1e8, 0, 1e8)
@time quadgk(f, -1e8, 0, 1e8)

# %%
dump(f)

# %%
a = 1.0
g = let a = a; x -> sin(a*x)/(a*x) end
@time quadgk(g, -1e8, 0, 1e8)
@time quadgk(g, -1e8, 0, 1e8)
@time quadgk(g, -1e8, 0, 1e8)

# %%
dump(g)

# %%
function F(a = 1.0)
    quadgk(x -> sin(a*x)/(a*x), -1e8, 0, 1e8)
end

@time F()
@time F()
@time F()

# %%
