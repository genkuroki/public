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
# Since there is no minimal working example (MWE), I could not figure out what you really wanted to do.
# (Please read https://discourse.julialang.org/t/please-read-make-it-easier-to-help-you/14757)
#
# However, I thought that possibly you would like to do something like this. Please take a look at the MWE below.

# %%
module SCI

export Problem, Algorithm_A, Algorithm_B, solve

abstract type AbstractProblem end
struct Problem{T} <: AbstractProblem x::T end

abstract type AbstractAlgorithm end
Base.@kwdef struct Algorithm_A{T} <: AbstractAlgorithm a::T = 2.0 end
Base.@kwdef struct Algorithm_B{T} <: AbstractAlgorithm a::T = 3.0 end
default_algorithm(prob::Problem) = Algorithm_A()

struct Solution{R, P<:AbstractProblem, A<:AbstractAlgorithm} result::R; prob::P; alg::A end
solve(prob::AbstractProblem) = solve(prob, default_algorithm(prob))
solve(prob::AbstractProblem, alg::AbstractAlgorithm) = Solution(alg.a * prob.x, prob, alg)

"""
Here, `Base` module plays the role of `COM` package.
SCI module can define the method of 
`Base.show` = `COM.show` function for `SCI.Solution` type.
"""
function Base.show(io::IO, sol::Solution)
    result = """
    Problem:   $(sol.prob)
    Algorithm: $(sol.alg)
    Result:    $(sol.result)
    """
    print(io, result)
end

end

using .SCI
@show prob = Problem(1.5)
println()
println(solve(prob))
println(solve(prob, Algorithm_B()))

# %%
module SCI2 # extension of SCI module

export Problem2, Algorithm_C

using ..SCI: SCI, AbstractProblem, AbstractAlgorithm, Solution

struct Problem2{T} <: AbstractProblem x::T; y::T end
Base.@kwdef struct Algorithm_C{T} <: AbstractAlgorithm a::T = 2.0; b = 3.0 end
SCI.default_algorithm(prob::Problem2) = Algorithm_C()
SCI.solve(prob::Problem2, alg::Algorithm_C) = Solution(alg.a * prob.x + alg.b * prob.y, prob, alg)

end

using .SCI2
@show prob2 = Problem2(10.0, 1.0)
println()
println(solve(prob2))

# %% [markdown]
# In the above example I don't define `AbstractSolution` type, nor did I make `Solution` a subtype of it. But you can define it, and then define another Solution type as a subtype of it in `SCI2` module, and also define how to display objects of another Solution type.   Many other extensions are also possible.
#
# See also https://discourse.julialang.org/t/function-depending-on-the-global-variable-inside-module/64322/10

# %%
