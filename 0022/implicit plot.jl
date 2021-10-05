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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using Plots
h(x, y) = x^(x^y) - y

x = range(0, 0.11; length=1000)
y = range(0, 1.1; length=1000)
contour(x, y, h; levels=zeros(1), color=:red, colorbar=false)
plot!(; xlim=extrema(x), ylim=extrema(y), size=(400, 400))

# %%
using Plots

function f(x, y, n)
    z = x^y
    for i in 1:n
        z = x^z
    end
    z
end

x = range(0, 2.1, length=1000)
y = range(0, 10.1; length=1000)
plot(; colorbar=false)
for n in 1:4
    contour!(x, y, (x, y) -> f(x, y, n) - y; levels=zeros(1), c=n, ls=:auto)
end
plot!(; xlim=extrema(x), ylim=extrema(y), size=(400, 400))

# %%
