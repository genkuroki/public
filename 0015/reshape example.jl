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
module Foo

const a = Float64[]

function prepare_workspace(m::Integer...)
    resize!(a, prod(m))
    reshape(a, m)
end

prepare_workspace(x::AbstractArray) = prepare_workspace(size(x)...)

function f(x)
    ws = prepare_workspace(x)
    ws .= 2x
end

end

@show Foo.f([1 2; 3 4])
@show Foo.a;

# %%
