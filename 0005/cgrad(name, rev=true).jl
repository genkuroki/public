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

# %%
using Plots
h=0.05
y, x = range(-5,  stop=5+h, step=h), range(-5, stop=10+h, step=h)
X = [xx for yy in y, xx in x]
Y = [yy for yy in y, xx in x]   
z = (sin.(X)) .^10 + cos.(10 .+ Y .* X) + cos.(X) + 0.2*Y + 0.1*X
heatmap(x, y, z, c=:viridis)

# %%
using Plots
h=0.05
y, x = range(-5,  stop=5+h, step=h), range(-5, stop=10+h, step=h)
X = [xx for yy in y, xx in x]
Y = [yy for yy in y, xx in x]   
z = (sin.(X)) .^10 + cos.(10 .+ Y .* X) + cos.(X) + 0.2*Y + 0.1*X
heatmap(x, y, z, c=reverse(cgrad(:viridis)))

# %%
using Plots
h=0.05
y, x = range(-5,  stop=5+h, step=h), range(-5, stop=10+h, step=h)
X = [xx for yy in y, xx in x]
Y = [yy for yy in y, xx in x]   
z = (sin.(X)) .^10 + cos.(10 .+ Y .* X) + cos.(X) + 0.2*Y + 0.1*X
heatmap(x, y, z, c=cgrad(:viridis, rev=true))

# %%
using Plots
h=0.05
x, y = -5:h:10+h, -5:h:5+h
f(x, y) = sin(x) ^10 + cos(10 + x * y) + cos(x) + 0.2y + 0.1x
heatmap(x, y, f, c=cgrad(:viridis, rev=true))

# %% [markdown]
# You can simply do
#
# ```julia
# heatmap(x, y, z, c=cgrad(:viridis, rev=true))
# ```
#
# or
#
# ```julia
# heatmap(x, y, z, c=reverse(cgrad(:viridis)))
# ```
#
# ref. [ColorGradient](http://docs.juliaplots.org/latest/colorschemes/#ColorGradient)
#
# P.S. More simply, 
#
# ```julia
# using Plots
# h=0.05
# x, y = -5:h:10+h, -5:h:5+h
# f(x, y) = sin(x) ^10 + cos(10 + x * y) + cos(x) + 0.2y + 0.1x
# heatmap(x, y, f, c=cgrad(:viridis, rev=true))
# ```
#
# A lot of dots in one line can often cause bugs

# %%
