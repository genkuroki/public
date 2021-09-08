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

# %% [markdown]
# # Attempt to use ConcreteStructs.jl with Parameters.jl
#
# * Gen Kuroki
# * 2021-09-06
# * https://github.com/jonniedie/ConcreteStructs.jl/issues/4#issuecomment-913998030
#
# I have tried to make ConcreteStructs.jl and Parameters.jl work together well.
#
# __Conclusion:__ It is possible to do so by making the following two changes.
#
# * Change `@concrete` not creating the inner constructor, so that the `Foo{__T_a, __T_b, __T_c}(a, b, c)`-type default constructor will be defined.
# * Change `@with_kw` expanding macros in the argument, [like `Base.@kwdef`](https://github.com/JuliaLang/julia/blob/4931faa34a8a1c98b39fb52ed4eb277729120128/base/util.jl#L455).
#
# Then `@concrete` works well with `@with_kw` and more completely with `Base.@kwdef`.
#
# See below for details.

# %%
VERSION

# %%
using ConcreteStructs
using Parameters

macro macroexpand_rmln(code)
    :(macroexpand($__module__, $(QuoteNode(code)), recursive=true) |>
        Base.remove_linenums!)
end

# %% [markdown]
# ## Plain struct

# %%
struct Foo{A, B, C}
    a::A
    b::B
    c::C
end

# %%
methods(Foo)

# %%
methods(Foo{1,2,3})

# %% [markdown]
# The default constructor `Foo{A, B, C}(a, b, c)` is defined.

# %% [markdown]
# ## @concrete struct
#
# `@concrete` removes the `Foo_concrete{__T_a, __T_b, __T_c}(a, b, c)`-type default constructor.

# %%
@concrete struct Foo_concrete
    a
    b
    c
end

# %%
methods(Foo_concrete)

# %%
methods(Foo_concrete{1,2,3})

# %% [markdown]
# `Foo_concrete{__T_a, __T_b, __T_c}(a, b, c)` is not defined because only the inner constructor `Foo_concrete(a::__T_a, b::__T_b, c::__T_c)` is defined.

# %%
@macroexpand_rmln @concrete struct Foo_concrete
    a
    b
    c
end

# %% [markdown]
# ## Base.@kwdef concrete struct

# %%
Base.@kwdef @concrete struct Foo_kwdef_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %%
Foo_kwdef_concrete()

# %%
Foo_kwdef_concrete{Int, Float64, String}(a = 4, b = 5.0, c = "six")

# %% [markdown]
# The reason for this error is that `Foo_kwdef_concrete{__T_a, __T_b, __T_c}(a, b, c)` is not defined.

# %%
methods(Foo_kwdef_concrete)

# %%
methods(Foo_kwdef_concrete{1,2,3})

# %%
@macroexpand_rmln Base.@kwdef @concrete struct Foo_kwdef_concrete
    a
    b
    c
end

# %% [markdown]
# ## @with_kw @concrete struct causes error

# %%
@with_kw @concrete struct Foo_with_kw_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %% [markdown]
# [The first line of `Base.@kwdef`](https://github.com/JuliaLang/julia/blob/4931faa34a8a1c98b39fb52ed4eb277729120128/base/util.jl#L455) expands macros in the argument expression:
#
# ```julia
#     expr = macroexpand(__module__, expr) # to expand @static
# ```
#
# This is what makes `Base.@kwdef @concrete struct` possible.  Make a similar change to `Parameters.@with_kw`.

# %% [markdown]
# ## Change @with_kw expanding macros in the argument

# %%
@eval Parameters macro with_kw(typedef)
    typedef = macroexpand(__module__, typedef) # inserted
    return esc(with_kw(typedef, __module__, true))
end

# %%
@with_kw @concrete struct Foo_with_kw_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %% [markdown]
# Okay, it seems to have worked.  But...

# %%
methods(Foo_with_kw_concrete)

# %%
methods(Foo_with_kw_concrete{1,2,3})

# %%
Foo_with_kw_concrete()

# %% [markdown]
# The reason for this error is that `Foo_with_kw_concrete{__T_a, __T_b, __T_c}(a, b, c)` is not defined again.

# %%
@macroexpand_rmln @with_kw @concrete struct Foo_with_kw_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %% [markdown]
# ## Change `@concrete` not creating the inner constructor
#
# Change `ConcreteStructs._concretize(expr)` not creating the inner constructor.

# %%
@eval ConcreteStructs function _concretize(expr)
    expr isa Expr && expr.head == :struct || error("Invalid usage of @concrete")
    
    is_mutable = expr.args[1]
    struct_name, type_params, super = _parse_head(expr.args[2])
    line_tuples = _parse_line.(expr.args[3].args)
    lines = first.(line_tuples)
    type_params_full = (type_params..., filter(x -> x!==nothing, last.(line_tuples))...)

    struct_type = if length(type_params_full) == 0
        struct_name
    else
        Expr(:curly, struct_name, type_params_full...)
    end

    head = Expr(:(<:), struct_type, super)
    # constructor_expr = _make_constructor(struct_name, type_params, type_params_full, lines)
    # body = Expr(:block, lines..., constructor_expr)
    body = Expr(:block, lines...)
    struct_expr = Expr(:struct, is_mutable, head, body)
    
    return struct_expr, struct_name, type_params
end

# %%
@macroexpand_rmln @concrete struct Bar_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %% [markdown]
# The inner constructor has been deleted.

# %% [markdown]
# ## @concrete works well with @with_kw

# %%
@with_kw @concrete struct Bar_with_kw_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %%
methods(Bar_with_kw_concrete)

# %%
methods(Bar_with_kw_concrete{1,2,3})

# %% [markdown]
# The default constructor `Bar_with_kw_concrete{__T_a, __T_b, __T_c}(a, b, c)` is implicitly defined.

# %%
Bar_with_kw_concrete()

# %%
Bar_with_kw_concrete(c = '3')

# %%
Bar_with_kw_concrete{Int, Float64, String}(a = 4, b = 5.0, c = "six")

# %%
Bar_with_kw_concrete(4, 5.0, "six")

# %% [markdown]
# `ConcreteStructs.@concrete` works well with `Parameters.@with_kw`!

# %%
@macroexpand_rmln @with_kw @concrete struct Bar_with_kw_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %% [markdown]
# ## @concrete works well more completely with Base.@kwdef

# %%
Base.@kwdef @concrete struct Bar_kwdef_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %%
Bar_kwdef_concrete()

# %%
Bar_kwdef_concrete(c = '3')

# %%
Bar_kwdef_concrete{Int, Float64, String}(a = 4, b = 5.0, c = "six") # not error

# %%
Bar_kwdef_concrete(4, 5.0, "six")

# %% [markdown]
# `ConcreteStructs.@concrete` works more completely well with `Base.@kwdef`!

# %%
@macroexpand_rmln Base.@kwdef @concrete struct Bar_kwdef_concrete
    a = 1
    b = 2.0
    c = "three"
end

# %%
