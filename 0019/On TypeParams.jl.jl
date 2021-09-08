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
# cf. https://github.com/genkuroki/public/blob/main/0019/ConcreteStructs.jl%20with%20Parameters.jl.ipynb

# %%
using TypeParams

macro macroexpand_rmln(code)
    :(macroexpand($__module__, $(QuoteNode(code)), recursive=true) |>
        Base.remove_linenums!)
end

# %%
@macroexpand_rmln @typeparams struct Foo_typedef
    a::{}
    b::{}
    c::{}
end

# %%
@macroexpand_rmln @typeparams struct Foo_typedef
    a::{}
    b::{}
    c::{<:AbstractString}
end

# %%
@macroexpand_rmln Base.@kwdef @typeparams struct Foo_typedef
    a::{} = 1
    b::{} = 2.0
    c::{<:AbstractString} = "three"
end

# %%
using Parameters

# %%
Base.@kwdef @typeparams struct Foo_typedef
    a::{} = 1
    b::{} = 2.0
    c::{<:AbstractString} = "three"
end

Foo_typedef()

# %%
@with_kw @typeparams struct Foo_with_jw
    a::{} = 1
    b::{} = 2.0
    c::{<:AbstractString} = "three"
end

# %% [markdown]
# cf. https://github.com/JuliaLang/julia/blob/4931faa34a8a1c98b39fb52ed4eb277729120128/base/util.jl#L455

# %%
@eval Parameters macro with_kw(typedef)
    typedef = macroexpand(__module__, typedef) # inserted
    return esc(with_kw(typedef, __module__, true))
end

# %%
@with_kw @typeparams struct Foo_with_kw
    a::{} = 1
    b::{} = 2.0
    c::{<:AbstractString} = "three"
end

Foo_with_kw()

# %%
@macroexpand_rmln @with_kw @typeparams struct Foo_with_kw
    a::{} = 1
    b::{} = 2.0
    c::{<:AbstractString} = "three"
end

# %%
