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

# %% [markdown]
# [quote="BambOoxX, post:12, topic:64568"]
# you define the `SCI` package such that it patches the `COM` package is this right ?
# [/quote]
#
# Yes.  It is a standard pattern.
#
# [quote="BambOoxX, post:12, topic:64568"]
# If so : What would happen if `COM` is not loaded prior to `SCI` .
# [/quote]
#
# Assume that, unlike the MWE I have shown above, both COM and SCI are packaged and can be loaded with `using COM` and `using SCI`.  Then the result of `using SCI; using COM` and the result of `using COM; using SCI` will be the same.
#
# In the following example, `COM` package has a function `prettify(x)` whose methods can be defined for user-defined types, and a function `prettyprint(x)` that uses `prettify(x)`.   `SCI` package makes use of them.
#
# __src/COM.jl of COM package:__
# ```julia
# """
# `COM.prettyprint(x)` prints nicely the object `x` for which `COM.prettify(x)` method is defined.
# """
# module COM
# prettify(x) = sprint(io -> show(io, "text/plain", x))
# prettyprint(io::IO, x) = print(io, prettify(x))
# prettyprint(x) = prettyprint(stdout, x)
# end
# ```
#
# __Example of COM__
# ```
# using COM
# COM.prettyprint([π, 2π])
# ```
# __Output:__
# ```
# 2-element Vector{Float64}:
#  3.141592653589793
#  6.283185307179586
# ```
#
# __src/SCI.jl of SCI package (with deps COM):__
# ```julia
# """A scirntific module (calculate exp(x))"""
# module SCI
#
# abstract type AbstractProblem end
# struct Problem{T} <: AbstractProblem x::T end
#
# abstract type AbstractAlgorithm end
# struct Builtin <: AbstractAlgorithm end
# Base.@kwdef struct Taylor <: AbstractAlgorithm n::Int = 10 end
# default_algorithm(prob::Problem) = Builtin()
#
# struct Solution{R, P<:AbstractProblem, A<:AbstractAlgorithm} result::R; prob::P; alg::A end
# solve(prob::AbstractProblem) = solve(prob, default_algorithm(prob))
# solve(prob::AbstractProblem, alg::Builtin) = Solution(exp(prob.x), prob, alg)
# solve(prob::AbstractProblem, alg::Taylor) = Solution(sum(prob.x^k/factorial(k) for k in 0:alg.n), prob, alg)
#
# using COM
#
# COM.prettify(sol::Solution{R, P, A}) where {R, P<:AbstractProblem, A<:Builtin} = """
# Problem:   x = $(sol.prob.x)
# Algorithm: builtin exp(x)
# Result:    $(sol.result)
# """
#
# COM.prettify(sol::Solution{R, P, A}) where {R, P<:AbstractProblem, A<:Taylor} = """
# Problem:   x = $(sol.prob.x)
# Algorithm: Taylor series of exp(x) upto degree $(sol.alg.n)
# Result:    $(sol.result)
# """
#
# end
# ```
#
# __Code to run:__
# ```julia
# # The order of the following two lines may be reversed.
# using SCI
# using COM
#
# prob = SCI.Problem(1)
# COM.prettyprint(SCI.solve(prob))
# println()
# COM.prettyprint(SCI.solve(prob, SCI.Taylor()))
# ```
#
# __Output:__
# ```
# Problem:   x = 1
# Algorithm: builtin exp(x)
# Result:    2.718281828459045
#
# Problem:   x = 1
# Algorithm: Taylor series of exp(x) upto degree 10
# Result:    2.7182818011463845
# ```
#
# In this way, you can use the pretty printing functions provided by `COM` package in any other package without changing any code in `COM` package.

# %%
# The order of the following two lines may be reversed.
using SCI
using COM

prob = SCI.Problem(1)
COM.prettyprint(SCI.solve(prob))
println()
COM.prettyprint(SCI.solve(prob, SCI.Taylor()))
