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
module O

struct Foo{A, B}
    a::A
    b::B
end

struct Bar{A, B}
    a::A
    b::B
    Bar(a::A = 12) where A = new{A, A}(a, 2a)
end

end

# %%
O.Foo(12, 34)

# %%
O.Foo{Float64, Float64}(12, 34)

# %%
methods(O.Foo)

# %%
methods(O.Foo{Float64, Float64})

# %%
O.Bar()

# %%
O.Bar(12.34)

# %%
methods(O.Bar)

# %%
methods(O.Bar{Float64, Float64})

# %%
O.Bar(Int8(3))

# %%
O.Bar(Int8(64))

# %%
@eval O struct Bar{A, B}
    a::A
    b::B
    Bar(a::A = 12) where A = new{A, A}(a, A(2)a)
end

# %%
O.Bar(Int8(64))

# %%
O.Bar(UInt8(128))

# %%
2UInt8(128)

# %%
UInt8(2)UInt8(128)

# %%
