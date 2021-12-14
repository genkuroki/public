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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %%
struct Foo{A, B}
    a::A
    b::B
end

# %%
methods(Foo)

# %%
methods(Foo{Int, Int})

# %%
foo = Foo(1.2, 34)

# %%
mult(x::Foo{<:Number, <:Number}) = x.a * x.b

# %%
methods(mult)

# %%
mult(Foo(2, 3 + 4im))

# %%
mult(x::Foo{<:Integer, <:AbstractString}) = x.b ^ x.a

# %%
methods(mult)

# %%
mult(Foo(10, "hoge! "))

# %%
