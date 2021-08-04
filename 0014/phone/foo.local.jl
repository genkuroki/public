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
struct Foo
    var"local"::Float64
end
@show foo = Foo(1.2)
@show foo.local;

# %%
:(struct Foo
    var"local"::Float64
end) |> Base.remove_linenums! |> print

# %%
