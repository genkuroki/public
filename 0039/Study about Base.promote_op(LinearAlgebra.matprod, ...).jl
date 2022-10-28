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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %%
using LinearAlgebra
using Polynomials
using SparseArrays

# %%
p = Polynomial([0//1,0//1,1//1], :h)

# %%
typeof(p)

# %%
p = Polynomial(Rational[0,0,1], :h)

# %%
a = sparse([p;;])

# %%
promote_type(eltype(a), eltype(-a))

# %%
a*a

# %%
-a

# %%
(-a)*a

# %%
@which (-a)*a

# %%
@which SparseArrays.spmatmul(-a, a)

# %%
Base.promote_op(LinearAlgebra.matprod, Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h})

# %%
@which Base.promote_op(LinearAlgebra.matprod, Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h})

# %%
Base._return_type(LinearAlgebra.matprod, Tuple{Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h}})

# %%
@which Base._return_type(LinearAlgebra.matprod, Tuple{Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h}})

# %%
Core.Compiler.return_type(LinearAlgebra.matprod, Tuple{Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h}})

# %%
@which Core.Compiler.return_type(LinearAlgebra.matprod, Tuple{Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h}})

# %%
@eval Core.Compiler function my_return_type(@nospecialize(f), t::DataType) # this method has a special tfunc
    @show world = ccall(:jl_get_tls_world_age, UInt, ())
    println()
    @show args = Any[_return_type, NativeInterpreter(world), Tuple{Core.Typeof(f), t.parameters...}]
    println()
    return @show ccall(:jl_call_in_typeinf_world, Any, (Ptr{Ptr{Cvoid}}, Cint), args, length(args))
end

# %%
Core.Compiler.my_return_type(LinearAlgebra.matprod, Tuple{Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h}})

# %%
Base.promote_op(LinearAlgebra.matprod, Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h})

# %%
