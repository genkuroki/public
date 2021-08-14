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
function f(L)
    a = []
    for i in 1:L
        push!(a, rand())
    end
    a
end

function g(L)
    a = Float64[]
    for i in 1:L
        push!(a, rand())
    end
    a
end

function h(L)
    a = Vector{Float64}(undef, L)
    for i in eachindex(a)
        a[i] = rand()
    end
    a
end

f(1); g(1); h(1); # compile

# %%
print("a = [] and push!:                           ")
@time a = f(10^6)
print("a = Float64[] and push!:                    ")
@time b = g(10^6)
print("a = Vector{Float64}(undef, L) and a[i] =...:")
@time c = h(10^6);

# %%
using BenchmarkTools

println("Test of constructors")
print("a = [] and push!:                          ")
@btime a = f(10^6)
print("a = Float64[] and push!:                    ")
@btime b = g(10^6)
print("a = Vector{Float64}(undef, L) and a[i] =...:")
@btime c = h(10^6);

# %%
println("Test of sum")
print("a = [] and push!:                        ")
@btime sum($a)
print("a = Float64[] and push!:                    ")
@btime sum($b)
print("a = Vector{Float64}(undef, L) and a[i] =...:")
@btime sum($c);

# %%
