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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
using BenchmarkTools
versioninfo()

# %%
function f(L=10^8)
    s = 0.0
    for k in 1:L
        s += sin(k)/k
    end
    1 + 2s
end

@btime f()

# %%
@code_native debuginfo=:none f()

# %%
using LoopVectorization

"single-thread optimized version"
function f_turbo(L=10^8)
    s = 0.0
    @turbo for k in 1:L
        s += sin(k)/k
    end
    1 + 2s
end

@btime f_turbo()

# %%
using LoopVectorization

"multi-thread optimized version"
function f_tturbo(L=10^8)
    s = 0.0
    @tturbo for k in 1:L
        s += sin(k)/k
    end
    1 + 2s
end

@show Threads.nthreads()
@btime f_tturbo()

# %%
@code_native debuginfo=:none f_turbo(10^8)

# %%
@macroexpand function f_turbo(L=10^8)
    s = 0.0
    @turbo for k in 1:L
        s += sin(k)/k
    end
    1 + 2s
end

# %%
@time f_tturbo(10^9)

# %%
using LoopVectorization

@inline g(x) = sin(x)/x

"single-thread optimization"
function f_turbo_inline(L=10^8)
    s = 0.0
    @turbo for k in 1:L
        s += g(k)
    end
    1 + 2s
end

@btime f_turbo_inline()

# %%
using LoopVectorization

@noinline g(x) = sin(x)/x

"single-thread optimization"
function f_turbo_noinline(L=10^8)
    s = 0.0
    @turbo for k in 1:L
        s += g(k)
    end
    1 + 2s
end

@btime f_turbo_noinline()

# %%
