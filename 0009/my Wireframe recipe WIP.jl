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
module O

using Plots

struct Wireframe{F, X, Y}
    f::F
    x::X
    y::Y
end

vecwithnan(X) = [
    vec([X ; fill(NaN, 1, size(X,  2))])
    vec([X'; fill(NaN, 1, size(X', 2))])
]

@recipe function f(wf::Wireframe)
    f, x, y = wf.f, wf.x, wf.y
    X, Y = reim(complex.(x', y))
    XX, YY = vecwithnan(X), vecwithnan(Y)
    ZZ = f.(XX, YY)
    XX, YY, ZZ
end

@recipe function f(wf::Wireframe{<:AbstractMatrix})
    Z, x, y = wf.f, wf.x, wf.y
    X, Y = reim(complex.(x', y))
    XX, YY, ZZ = vecwithnan(X), vecwithnan(Y), vecwithnan(Z)
    XX, YY, ZZ
end

end

# %%
using Plots
x = y = range(-1, 1, 51)
f(x, y) = x == y == 0 ? zero(x) : x^2*y/(x^2 + y^2)
plot(O.Wireframe(f, x, y); label="")

# %%
Z = f.(x', y)
plot(O.Wireframe(Z, x, y); label="")

# %%
wireframe(x, y, f)

# %%
