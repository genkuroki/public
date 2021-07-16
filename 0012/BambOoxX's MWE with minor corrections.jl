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

# %% [markdown]
# https://discourse.julialang.org/t/handling-of-package-inter-dependency/64568/24

# %%
module SCI
export fun1,fun2
fun1(x) = x
fun2(x) = x^2
end

module COM
using ..SCI
"""function to initialize a container"""
init_container() = []
"""function to read a message from a pipe, assume it reads from the argument for simplicity"""
read_message(x) = x
# function to check if instruction is authorized (omitted for simplicity)

"""function doing something based on the message received"""
function compute!(container, msg)
    if msg == :fun1
        push!(container, fun1(rand()))
    elseif msg ==:fun2
        push!(container, fun2(rand()))
    end
    return container
end
end

# Load COM+SCI utilities
using .COM

@show container = COM.init_container()
@show msg = COM.read_message(:fun1) 
@show container = COM.compute!(container, msg)
@show msg = COM.read_message(:fun2) 
@show container = COM.compute!(container, msg);

# %%
