# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %%
?propertynames

# %%
module O

Base.@kwdef mutable struct Foo
    pub::String = "meow"
    priv::Int = 3
end

function Base.getproperty(foo::Foo, name::Symbol)
    name === :multiple && return getfield(foo, :pub)^getfield(foo, :priv)
    hasproperty(foo, name) && return getfield(foo, name)
    error("type Foo has no public property $name")
end

Base.show(io::IO, foo::Foo) = print(io, "Foo(pub = ", repr(foo.pub), ')')

Base.propertynames(foo::Foo, private::Bool=false) =
    private ? (fieldnames(typeof(foo))..., :multiple) : (:pub, :multiple)

end

# %%
foo = O.Foo()

# %%
propertynames(foo)

# %%
@which propertynames(foo)

# %%
methodswith(O.Foo)

# %%
propertynames(foo, true)

# %%
foo.multiple

# %%
foo.pub = "bow"
foo

# %%
foo.multiple

# %%
foo.priv

# %%
foo.pub

# %%
