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
# https://nbviewer.jupyter.org/urls/gist.githubusercontent.com/terasakisatoshi/e4af4f618161d87cdc288285585527c5/raw/3a50de35748be12ab6461156440685052aadfcd2/hexagon_something.ipynb

# %%
using HTTP, JSON
function display_tweet(link)
    api = "https://publish.twitter.com/oembed?url=$link"
    r = response = HTTP.request("GET", api);
    j = JSON.parse(String(r.body))
    HTML(j["html"])
end

# %% [markdown]
# https://twitter.com/dannchu/status/1443524638252810248

# %%
display_tweet("https://twitter.com/dannchu/status/1443524638252810248")

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
QQ = []
for n in 1:10
    @time countf!(X, a, n, L); flush(stdout)
    Q = bar(1:21, X/L; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(QQ, Q)
end
plot(QQ...; layout=(5, 2), size=(700, 1000))

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
QQ = []
for n in 1:10
    @time g!(X, a, n); flush(stdout)
    Q = bar(1:21, X/6^n; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(QQ, Q)
end
plot(QQ...; layout=(5, 2), size=(700, 1000))

# %%
X = zeros(Int, 21)
a = falses(6)
QQ = []
for n in 1:10
    @time g!(X, a, n); flush(stdout)
    Q = bar(1:21, X/6^n; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(QQ, Q)
end
plot(QQ...; layout=(5, 2), size=(700, 1000))

# %%
X = zeros(Int, 21)
a = falses(6)
QQ = []
for n in 1:10
    iter = Iterators.product(ntuple(_ -> 1:6, n)...)
    @time g!(X, a, iter); flush(stdout)
    Q = bar(1:21, X/6^n; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(QQ, Q)
end
plot(QQ...; layout=(5, 2), size=(700, 1000))

# %%
@code_warntype g!(X, a, Iterators.product(ntuple(_ -> 1:6, 10)...))

# %% [markdown]
# https://twitter.com/croce1/status/1445311359973937155

# %%
display_tweet("https://twitter.com/croce1/status/1445311359973937155")

# %%
using Combinatorics
using Memoization

P(n, s) = sum(k -> p(n, k) * q(n, k, s), 1:min(6, n))
@memoize p(n, k) = c(n, k)//6^n
c(n, k) = binomial(6, k) * sum(i -> binomial(k, i)*(-1)^(i+k)*i^n, 1:k)
S(k) = (a = binomial(k+1, 2); a:(7k - a))
@memoize B(k, s) = Iterators.filter(A -> sum(A) == s, powerset(1:6, k, k))
@memoize b(k, s) = count(_ -> true, B(k, s))
@memoize q(n, k, s) = b(k, s)//sum(t -> b(k, t), S(k))

X = zeros(21)
QQ = []
for n in 1:10
    @time for s in 1:21
        X[s] = P(n, s)
    end
    Q = bar(1:21, X; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(QQ, Q)
end
plot(QQ...; layout=(5, 2), size=(700, 1000))

# %%
X = zeros(21)
QQ = []
for n in 1:10
    @time for s in 1:21
        X[s] = P(n, s)
    end
    Q = bar(1:21, X; alpha=0.3, xtick=1:21, label="", title="n = $n")
    push!(QQ, Q)
end
plot(QQ...; layout=(5, 2), size=(700, 1000))

# %%
@memoize G(n) = g!(zeros(Int, 21), falses(6), n)
g(n, s) = G(n)[s]//6^n

ENV["COLUMNS"] = 130
@time M = [g(n, s) - P(n, s) for n in 1:10, s in 1:21]

# %%
