# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %%
using BenchmarkTools

# %%
using Random

function mcpi(N, rng = MersenneTwister())
    c = 0
    for _ in 1:N
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/N
end

println("Julia v", VERSION)
print("mcpi(10^8):")
@btime mcpi(10^8)

# %%
using Base.Threads
using Random

function mcpi_atomic(N, rng = MersenneTwister())
    a = Atomic{Int}(0)
    @threads for _ in 1:N
        atomic_add!(a, Int(rand(rng)^2 + rand(rng)^2 ≤ 1))
    end
    4a[]/N
end

println("Julia v", VERSION)
@show Threads.nthreads()
print("mcpi_atomic(10^8):")
@btime mcpi_atomic(10^8)

# %% [markdown]
# https://gist.github.com/TsuMakoto/4138d3c2fd05a400d02eff0b91e34658

# %%
using Base.Threads
using Random

mutable struct MCAtomic{T}; @atomic n::T; end

function mcpi_tm(N, rng = MersenneTwister())
    mc = MCAtomic(0)
    @threads for _ in 1:N
        @atomic mc.n += Int(rand(rng)^2 + rand(rng)^2 ≤ 1)
    end
    4mc.n/N
end

println("Julia v", VERSION)
@show Threads.nthreads()
print("mcpi_tm(10^8):")
@btime mcpi_tm(10^8)

# %%
using Base.Threads
using Distributed: splitrange
using Random

function mcpi_splittange_atomic(N)
    ranges = splitrange(1, N, nthreads())
    a = Atomic{Int}(0)
    @threads for ran in ranges
        rng = MersenneTwister()
        c = 0
        for _ in ran
            c += rand(rng)^2 + rand(rng)^2 ≤ 1
        end
        atomic_add!(a, c)
    end
    4a[]/N
end

println("Julia v", VERSION)
@show Threads.nthreads()
print("mcpi_splittange_atomic(10^8):")
@btime mcpi_splittange_atomic(10^8)

# %%
# The following code is a modified version of
# 
#    function _threadsfor(iter, lbody, schedule) 
#    macro threads(args...)
#
# in https://github.com/JuliaLang/julia/blob/9f3265399227fbfc4f0160ec3592a9262bd3eb5f/base/threadingconstructs.jl
#
# Its license is MIT: https://julialang.org/license

using Base.Threads
using Base.Threads: threading_run

function _my_threadsfor(iter, lbody, prebody, postbody, schedule)
    lidx = iter.args[1]         # index
    range = iter.args[2]
    quote
        local threadsfor_fun
        let range = $(esc(range))
        function threadsfor_fun(onethread=false)
            r = range # Load into local variable
            lenr = length(r)
            # divide loop iterations among threads
            if onethread
                tid = 1
                len, rem = lenr, 0
            else
                tid = threadid()
                len, rem = divrem(lenr, nthreads())
            end
            # not enough iterations for all the threads?
            if len == 0
                if tid > rem
                    return
                end
                len, rem = 1, 0
            end
            # compute this thread's iterations
            f = firstindex(r) + ((tid-1) * len)
            l = f + len - 1
            # distribute remaining iterations evenly
            if rem > 0
                if tid <= rem
                    f = f + (tid-1)
                    l = l + tid
                else
                    f = f + rem
                    l = l + rem
                end
            end
            # run prebody
            $(esc(prebody))
            # run this thread's iterations
            for i = f:l
                local $(esc(lidx)) = @inbounds r[i]
                $(esc(lbody))
            end
            # run postbody
            $(esc(postbody))
        end
        end
        if threadid() != 1 || ccall(:jl_in_threaded_region, Cint, ()) != 0
            $(if schedule === :static
              :(error("`@my_threads :static` can only be used from thread 1 and not nested"))
              else
              # only use threads when called from thread 1, outside @threads
              :(Base.invokelatest(threadsfor_fun, true))
              end)
        else
            threading_run(threadsfor_fun)
        end
        nothing
    end
end

"""
    @my_threads
A macro to parallelize a `for` loop to run with multiple threads. 
It splits the iteration space among multiple tasks with `prebody` and `postbody`.
It runs those tasks on threads according to a scheduling policy.
Usage:
```julia
@my_threads [schedule] begin
    prebody
end for ...
    ...
end begin
    postbody
end
```
"""
macro my_threads(args...)
    na = length(args)
    if na == 4
        sched, prebody, ex, bostbody = args
        if sched isa QuoteNode
            sched = sched.value
        elseif sched isa Symbol
            # for now only allow quoted symbols
            sched = nothing
        end
        if sched !== :static
            throw(ArgumentError("unsupported schedule argument in @threads"))
        end
    elseif na == 3
        sched = :default
        prebody, ex, postbody = args
    else
        throw(ArgumentError("wrong number of arguments in @my_threads"))
    end
    if !(isa(ex, Expr) && ex.head === :for)
        throw(ArgumentError("@my_threads requires a `for` loop expression"))
    end
    if !(ex.args[1] isa Expr && ex.args[1].head === :(=))
        throw(ArgumentError("nested outer loops are not currently supported by @my_threads"))
    end
    return _my_threadsfor(ex.args[1], ex.args[2], prebody, postbody, sched)
end

@doc @my_threads

# %%
using Random

function mcpi_my_threads(N)
    a = Atomic{Int}(0)
    @my_threads begin
        rng = MersenneTwister()
        c = 0
    end for _ in 1:N
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end begin
        atomic_add!(a, c)
    end
    4a[]/N
end

println("Julia v", VERSION)
@show Threads.nthreads()
print("mcpi_my_threads(10^8):")
@btime mcpi_my_threads(10^8)

# %%
