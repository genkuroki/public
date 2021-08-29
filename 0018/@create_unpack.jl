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
"""
The `CreateUnPack` module only exports the `@create_unpack` macro,
which creates a macro unpacking fields or properties of an object.

__Simple Intended Usage:__
```
julia> using .CreateUnPack

julia> struct Foo{A, B, C} a::A; b::B; c::C end

julia> @create_unpack all_Foo Foo
@unpackall_Foo (macro with 1 method)

help?> @unpackall_Foo
@unpackall_Foo(x) unpacks the fields (:a, :b, :c) of x.

julia> @macroexpand @unpackall_Foo foo
:((a, b, c) = (foo.a, foo.b, foo.c))

julia> @unpackall_Foo Foo(1, 2.0, "three")
(1, 2.0, "three")

julia> a, b, c
(1, 2.0, "three")
```
"""
module CreateUnPack

export @create_unpack

"""
    @create_unpack(name::Symbol, expr)

creates the macro named `Symbol(:unpack, name)`
which unpacks the fields specified by `expr` of an object.

Let `val` be the value of `expr`.

* If `val` is a tuple of symbols, then set the unpacking fields to `val`.
* If `val` is a typp, then set the unpacking fields to `fieldnames(val)`.
* Otherwise, set the unpacking fields to `propertynames(val)`.

__Example:__

```
julia> @create_unpack _cat_and_dog (:cat, :dog)
@unpack_cat_and_dog (macro with 1 method)

help?> @unpack_cat_and_dog
@unpack_cat_and_dog(obj) unpacks the fields (:cat, :dog) of obj.

julia> @macroexpand @unpack_cat_and_dog x
:((cat, dog) = (x.cat, x.dog))

julia> @unpack_cat_and_dog (dog = "bowwow", mouse = "squeak", cat = "meow")
("meow", "bowwow")

julia> cat, dog
("meow", "bowwow")

julia> struct Foo{A, B, C} a::A; b::B; c::C end

julia> @create_unpack all_Foo Foo
@unpackall_Foo (macro with 1 method)

help?> @unpackall_Foo
@unpackall_Foo(obj) unpacks the fields (:a, :b, :c) of obj.

julia> @macroexpand @unpackall_Foo foo
:((a, b, c) = (foo.a, foo.b, foo.c))

julia> @unpackall_Foo Foo(1, 2.0, "three")
(1, 2.0, "three")

julia> a, b, c
(1, 2.0, "three")

julia> baz = (p = "one", q = 2.0, r = 3)

julia> @create_unpack all_baz baz
@unpackall_baz (macro with 1 method)

help?> @unpackall_baz
@unpackall_baz(obj) unpacks the fields (:p, :q, :r) of obj.

julia> @macroexpand @unpackall_baz baz
:((p, q, r) = (baz.p, baz.q, baz.r))

julia> @unpackall_baz baz
("one", '2', 3)

julia> p, q, r
("one", '2', 3)

```
"""
macro create_unpack(name::Symbol, expr)
    macroname = Symbol(:unpack, name)
    atmacroname = Symbol('@', macroname)
    val = Core.eval(__module__, expr)
    names = val isa NTuple{N, Symbol} where N ? val :
            val isa Type ? fieldnames(val) : propertynames(val)
    docstr = """`$atmacroname(obj)` unpacks the fields `$names` of `obj`."""
    quote
        macro $macroname(obj)
            Expr(:(=),
                Expr(:tuple, $names...),
                Expr(:tuple, (:($obj.$name) for name in $names)...)
            ) |> esc
        end
        @doc $docstr $(:($atmacroname))
        $atmacroname
    end |> esc
end

end

# %%
@doc CreateUnPack

# %%
using .CreateUnPack
@doc @create_unpack

# %%
@create_unpack _cat_and_dog (:cat, :dog)

# %%
?@unpack_cat_and_dog

# %%
@macroexpand @unpack_cat_and_dog x

# %%
@unpack_cat_and_dog (dog = "bowwow", mouse = "squeak", cat = "meow")

# %%
cat, dog

# %% tags=[]
struct Foo{A, B, C} a::A; b::B; c::C end

# %%
@create_unpack all_Foo Foo

# %%
?@unpackall_Foo

# %%
@macroexpand @unpackall_Foo foo

# %%
@unpackall_Foo Foo(1, 2.0, "three")

# %%
a, b, c

# %%
function f(foo::Foo)
    @unpackall_Foo foo
    @show a b c
    a, b, c
end

# %%
f(Foo(1, 2.0, "three"))

# %%
baz = (p = "one", q = 2.0, r = 3)

# %%
@create_unpack all_baz baz

# %%
?@unpackall_baz

# %%
@macroexpand @unpackall_baz baz

# %%
@unpackall_baz baz

# %%
p, q, r

# %%
module A
using ..CreateUnPack
struct Bar{X, Y} x::X; y::Y end
@create_unpack Bar Bar
end

# %%
?A.@unpackBar

# %%
@macroexpand A.@unpackBar bar

# %%
A.@unpackBar A.Bar(1, 2)

# %%
x, y

# %%
