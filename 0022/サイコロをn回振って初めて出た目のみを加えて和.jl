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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://twitter.com/dannchu/status/1443524638252810248

# %%
using Plots
default(titlefontsize=12, tickfontsize=6)

function f!(a, n)
    a .= false
    s = 0
    @inbounds for _ in 1:n
        i = rand(1:6)
        a[i] && continue
        a[i] = true
        s += i
        all(a) && break
    end
    s
end

function countf!(X, a, n, L)
    X .= 0
    for _ in 1:L
        @inbounds X[f!(a, n)] +=1
    end
    X
end

X = zeros(Int, 21)
a = falses(6)
L = 10^6
PP = []
for n in 1:10
    @time countf!(X, a, n, L); flush(stdout)
    P = bar(1:21, X/L; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(PP, P)
end
plot(PP...; layout=(5, 2), size=(700, 1000))

# %%
using Plots
default(titlefontsize=12, tickfontsize=6)

function g!(X, a, iter::Iterators.ProductIterator)
    X .= 0
    @inbounds for t in iter
        a .= false
        s = 0
        for i in t
            a[i] && continue
            a[i] = true
            s += i
            all(a) && break
        end
        X[s] += 1
    end
    X
end

function g!(X, a, n)
    iter = Iterators.product(ntuple(_ -> 1:6, n)...)
    g!(X, a, iter)
end

X = zeros(Int, 21)
a = falses(6)
PP = []
for n in 1:10
    @time g!(X, a, n); flush(stdout)
    P = bar(1:21, X/6^n; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(PP, P)
end
plot(PP...; layout=(5, 2), size=(700, 1000))

# %%
X = zeros(Int, 21)
a = falses(6)
PP = []
for n in 1:10
    @time g!(X, a, n); flush(stdout)
    P = bar(1:21, X/6^n; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(PP, P)
end
plot(PP...; layout=(5, 2), size=(700, 1000))

# %%
X = zeros(Int, 21)
a = falses(6)
PP = []
for n in 1:10
    iter = Iterators.product(ntuple(_ -> 1:6, n)...)
    @time g!(X, a, iter); flush(stdout)
    P = bar(1:21, X/6^n; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(PP, P)
end
plot(PP...; layout=(5, 2), size=(700, 1000))

# %%
@code_warntype g!(X, a, Iterators.product(ntuple(_ -> 1:6, 10)...))

# %%
