"""A scirntific module (calculate exp(x))"""
module SCI

abstract type AbstractProblem end
struct Problem{T} <: AbstractProblem x::T end

abstract type AbstractAlgorithm end
struct Builtin <: AbstractAlgorithm end
Base.@kwdef struct Taylor <: AbstractAlgorithm n::Int = 10 end
default_algorithm(prob::Problem) = Builtin()

struct Solution{R, P<:AbstractProblem, A<:AbstractAlgorithm} result::R; prob::P; alg::A end
solve(prob::AbstractProblem) = solve(prob, default_algorithm(prob))
solve(prob::AbstractProblem, alg::Builtin) = Solution(exp(prob.x), prob, alg)
solve(prob::AbstractProblem, alg::Taylor) = Solution(sum(prob.x^k/factorial(k) for k in 0:alg.n), prob, alg)

using COM

COM.prettify(sol::Solution{R, P, A}) where {R, P<:AbstractProblem, A<:Builtin} = """
Problem:   x = $(sol.prob.x)
Algorithm: builtin exp(x)
Result:    $(sol.result)
"""

COM.prettify(sol::Solution{R, P, A}) where {R, P<:AbstractProblem, A<:Taylor} = """
Problem:   x = $(sol.prob.x)
Algorithm: Taylor series of exp(x) upto degree $(sol.alg.n)
Result:    $(sol.result)
"""

end
