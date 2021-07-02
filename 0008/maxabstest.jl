# maxabstest.jl

using BenchmarkTools

N = isempty(ARGS) ? 500 : parse(Int, ARGS[1])

function max_abs(a, b)
    m = zero(promote_type(eltype(a), eltype(b)))
    for i in eachindex(a, b)
        tmp = abs(a[i] - b[i]) 
        tmp > m && (m = tmp)
    end 
    m
end

a = randn(N, N)
b = randn(N, N)

@show VERSION
@show Threads.nthreads()
println()
@show N
print("simple for loop:     ")
@btime max_abs($a, $b)
print("mapreduce abs∘- max: ")
@btime mapreduce(abs∘-, max, $a, $b)
print("maximum(abs, a - b): ")
@btime maximum(abs, $a - $b)
print("maximum(generator):  ")
@btime maximum(abs(i - j) for (i, j) in zip($a, $b))
print("maximum splat(abs∘-):")
@btime maximum(Base.splat(abs∘-), zip($a, $b))
using Tullio
print("Tullio:              ")
@btime @tullio (max) _ := abs($a[i] - $b[i])
print("Tullio (LoopVect.):  ")
using LoopVectorization
@btime @tullio (max) _ := abs($a[i] - $b[i])

function max_abs_turbo(a, b)
    m = zero(promote_type(eltype(a), eltype(b)))
    @turbo for i in eachindex(a, b)
        m = max(m, abs(a[i] - b[i]))
    end 
    m
end

function max_abs_tturbo(a, b)
    m = zero(promote_type(eltype(a), eltype(b)))
    @tturbo for i in eachindex(a, b)
        m = max(m, abs(a[i] - b[i]))
    end 
    m
end

print("LoopVect. @turbo:    ")
@btime max_abs_turbo($a, $b)
print("LoopVect. @tturbo:   ")
@btime max_abs_tturbo($a, $b)
