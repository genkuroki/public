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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
module COM

abstract type AbstractContainer end
struct Container{T<:Vector} <: AbstractContainer vec::T end
Base.parent(x::AbstractContainer) = getfield(x, :vec)
Base.pop!(x::AbstractContainer) = pop!(parent(x))
Base.push!(x::AbstractContainer, a...) = push!(parent(x), a...)

abstract type AbstractMessage end
struct Message{T} <: AbstractMessage instruction::T end
instruction(x::AbstractMessage) = getfield(x, :instruction)

init_container() = Container([])
read_message(x) = Message(x)

"""
function doing something based on the message received
compute! function will be defined elsewhere
"""
process!(compute!, container::AbstractContainer, msg::AbstractMessage) =
    compute!(container, msg)
end

module SCI

fun1(x) = x
fun2(x) = x^2

using ..COM

function compute!(container, msg)
    if msg == :fun1
        push!(container, fun1(rand()))
    elseif msg ==:fun2
        push!(container, fun2(rand()))
    end
    return container
end
end

@show compute! = SCI.compute!
@show container = COM.init_container()
@show msg = COM.read_message(:fun1) 
@show COM.process!(compute!, container, msg)
@show msg = COM.read_message(:fun2) 
@show COM.process!(compute!, container, msg);

# %%
module O
main(mod::Module) = print(mod.main())
end

module P
main() = float(π)
end

O.main(P);

# %%
