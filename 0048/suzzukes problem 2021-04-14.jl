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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# https://x.com/dannchu/status/1779438301398892932

# %%
using Distributions

function prob_orig(n)
    ab = rand(Bernoulli(0.6), n)
    ac = rand(Bernoulli(0.75), n)
    bc = rand(Bernoulli(0.6), n)
    k = 0
    for i in 1:n
        if [ab[i], ac[i], bc[i]] == [1, 0, 1] || [ab[i], ac[i], bc[i]] == [0, 1, 0]
            k += 1
        end
    end
    k/n
end

@time prob_orig(10^8)
@time prob_orig(10^8)
@time prob_orig(10^8)

# %%
using Distributions

function prob_rev1(n)
    k = 0
    for i in 1:n
        ab = rand(Bernoulli(0.6))
        ac = rand(Bernoulli(0.75))
        bc = rand(Bernoulli(0.6))
        if (ab, ac, bc) == (1, 0, 1) || (ab, ac, bc) == (0, 1, 0)
            k += 1
        end
    end
    k/n
end

@time prob_rev1(10^8)
@time prob_rev1(10^8)
@time prob_rev1(10^8)

# %%
using Distributions

function prob_rev2(n)
    k = 0
    for i in 1:n
        ab = rand(Bernoulli(0.6))
        ac = rand(Bernoulli(0.75))
        bc = rand(Bernoulli(0.6))
        k += ab == bc !== ac
    end
    k/n
end

@time prob_rev2(10^8)
@time prob_rev2(10^8)
@time prob_rev2(10^8)

# %%
using Distributions

function prob_rev3(n)
    k = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:n
        tid = Threads.threadid()
        ab = rand(Bernoulli(0.6))
        ac = rand(Bernoulli(0.75))
        bc = rand(Bernoulli(0.6))
        k[tid] += ab == bc !== ac
    end
    sum(k)/n
end

@time prob_rev3(10^8)
@time prob_rev3(10^8)
@time prob_rev3(10^8)

# %%
