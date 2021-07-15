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
using Nemo
matrix(T, A) = MatrixSpace(T, size(A)...)(A)
pl(x...) = println(x...)
macro s(expr)
    :(println($(string(expr)), " ="); show(stdout, "text/plain", $expr); println("\n"))
end

F, a = FiniteField(7, 2, "a")
P, x = F["x"]
X = matrix(F, (0:6)' .* a .+ (0:6))
A = F[a+1 a+2; a+3 3a+4]
v = F[5a+2; 3a+6]

@show F P; pl()
@show [F(3)^k for k in 1:6]; pl()
@show [a^k for k in 1:48]; pl()
@s X
@s rref(X)[2]
@show det(X) tr(X) rank(X) charpoly(P, X) minpoly(P, X); pl()
@s A
@s A'
@s A'A
@show det(A) tr(A) rank(A) charpoly(P, A) lu(A) inv(A); pl()
@show v A*v; pl()

# %%
methodswith(typeof(A))

# %%
