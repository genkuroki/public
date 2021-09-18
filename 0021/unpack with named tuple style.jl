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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %%
function f(x)
    (; a, b, c, d) = x
    @show a b c d
end

# %%
a = 1
b = 2.0
c = "three"
d = 0x04

nt = (; a, b, c, d)
@show nt;

# %%
f(nt);

# %%
function f(x)
    (; a, b, c, d) = x
    @show a b c d
end

# %%
Base.@kwdef struct Foo{A, B, C, D}
    a::A = 1
    b::B = 2.0
    c::C = "three"
    d::D = 0x04
end

foo = Foo()
@show foo;

# %%
f(foo);

# %%
function f(x)
    (; a, b, c, d) = x
    @show a b c d
end

# %%
using DataFrames

a = collect(1:5)
b = collect(1.0:5.0)
c = ["one", "two", "three", "four", "five"]
d = collect(0x01:0x05)

df = DataFrame(; a, b, c, d)

# %%
f(df);

# %%
@show nt;

# %%
typeof(nt)

# %%
dump(nt)

# %%
@show foo;

# %%
dump(foo)

# %%
g(x::NamedTuple{(:a, :b, :c, :d)}) = @show x
g(x::NamedTuple{(:a, :b, :c)}) = @show x

# %%
g(nt)

# %%
@which g(nt)

# %%
nt2 = (; a = 10, b = 20.0, c = "thirty")

# %%
g(nt2)

# %%
@which g(nt2)

# %%
