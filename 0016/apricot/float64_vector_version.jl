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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# Modified version of https://github.com/genkuroki/public/blob/main/0016/apricot/julia_translation_of_python_reimpl.ipynb

# %%
#using Seaborn
using ScikitLearn: @sk_import
@sk_import datasets: fetch_covtype
#using Random
#using StatsBase: sample

# %%
digits_data = fetch_covtype()

# %%
X_digits = permutedims(abs.(digits_data["data"]))
summary(X_digits)

# %%
"""`calculate_gains!(X, gains, current_values, idxs, current_concave_values_sum)` mutates `gains` only"""
function calculate_gains!(X, gains, current_values, idxs, current_concave_values_sum)
    Threads.@threads for i in eachindex(idxs)
        @inbounds idx = idxs[i]
        @inbounds gains[i] = sum(j -> sqrt(current_values[j] + X[j, idx]), axes(X, 1))
    end
    gains .-= current_concave_values_sum
end

@doc calculate_gains!

# %%
"""
The revised version of using popat!(idxs, idx) without mask
"""
function fit_popat(X, k; calculate_gains! = calculate_gains!)
    d, n = size(X)

    cost = 0.0

    ranking = Int[]
    total_gains = Float64[]

    current_values = zeros(d)
    current_concave_values_sum = sum(sqrt, current_values)

    idxs = collect(1:n)

    gains = zeros(n)
    while cost < k
        calculate_gains!(X, gains, current_values, idxs, current_concave_values_sum)

        idx = argmax(gains)
        best_idx = idxs[idx]
        curr_cost = 1.0
        
        cost + curr_cost > k && break

        cost += curr_cost
        # Calculate gains
        gain = gains[idx] * curr_cost

        # Select next
        current_values .+= @view X[:, best_idx]
        current_concave_values_sum = sum(sqrt, current_values)

        push!(ranking, best_idx)
        push!(total_gains, gain)

        popat!(idxs, idx)
    end
    return ranking, total_gains
end

# %%
"""
The revised version of using mask::BitVector with findall(mask)
"""
function fit_bitvector(X, k; calculate_gains! = calculate_gains!)
    d, n = size(X)

    cost = 0.0

    ranking = Int[]
    total_gains = Float64[]

    mask = trues(n) # `false` stands for "masked".
    current_values = zeros(d)
    current_concave_values_sum = sum(sqrt, current_values)

    idxs = collect(1:n)

    gains = zeros(n)
    while cost < k
        calculate_gains!(X, gains, current_values, idxs, current_concave_values_sum)

        idx = argmax(gains)
        best_idx = idxs[idx]
        curr_cost = 1.0
        
        cost + curr_cost > k && break

        cost += curr_cost
        # Calculate gains
        gain = gains[idx] * curr_cost

        # Select next
        current_values .+= @view X[:, best_idx]
        current_concave_values_sum = sum(sqrt, current_values)

        push!(ranking, best_idx)
        push!(total_gains, gain)

        mask[best_idx] = 0
        idxs = findall(mask)
    end
    return ranking, total_gains
end

# %%
"""
The revised version of using mask::Vector{Float64} with findall(mask .== 0)
"""
function fit_f64vector(X, k; calculate_gains! = calculate_gains!)
    d, n = size(X)

    cost = 0.0

    ranking = Int[]
    total_gains = Float64[]

    mask = zeros(n)
    current_values = zeros(d)
    current_concave_values_sum = sum(sqrt, current_values)

    idxs = collect(1:n)

    gains = zeros(n)
    while cost < k
        calculate_gains!(X, gains, current_values, idxs, current_concave_values_sum)

        idx = argmax(gains)
        best_idx = idxs[idx]
        curr_cost = 1.0
        
        cost + curr_cost > k && break

        cost += curr_cost
        # Calculate gains
        gain = gains[idx] * curr_cost

        # Select next
        current_values .+= @view X[:, best_idx]
        current_concave_values_sum = sum(sqrt, current_values)

        push!(ranking, best_idx)
        push!(total_gains, gain)

        mask[best_idx] = 1
        idxs = findall(mask .== 0)
    end
    return ranking, total_gains
end

# %%
"""
The original version of using mask::Vector{Float64} with findall(==(0), mask)
"""
function fit_f64vector_org(X, k; calculate_gains! = calculate_gains!)
    d, n = size(X)

    cost = 0.0

    ranking = Int[]
    total_gains = Float64[]

    mask = zeros(n)
    current_values = zeros(d)
    current_concave_values_sum = sum(sqrt, current_values)

    idxs = collect(1:n)

    gains = zeros(n)
    while cost < k
        calculate_gains!(X, gains, current_values, idxs, current_concave_values_sum)

        idx = argmax(gains)
        best_idx = idxs[idx]
        curr_cost = 1.0
        
        cost + curr_cost > k && break

        cost += curr_cost
        # Calculate gains
        gain = gains[idx] * curr_cost

        # Select next
        current_values .+= @view X[:, best_idx]
        current_concave_values_sum = sum(sqrt, current_values)

        push!(ranking, best_idx)
        push!(total_gains, gain)

        mask[best_idx] = 1
        idxs = findall(==(0), mask)
    end
    return ranking, total_gains
end

# %%
k = 1000

# %%
@time ranking0_pa, gains0_pa = fit_popat(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_pa, gains0_pa = fit_popat(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_pa, gains0_pa = fit_popat(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_bv, gains0_bv = fit_bitvector(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_bv, gains0_bv = fit_bitvector(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_bv, gains0_bv = fit_bitvector(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_f64v, gains0_f64v = fit_f64vector(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_f64v, gains0_f64v = fit_f64vector(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_f64v, gains0_f64v = fit_f64vector(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_org, gains0_org = fit_f64vector_org(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_org, gains0_org = fit_f64vector_org(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0_org, gains0_org = fit_f64vector_org(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@show ranking0_pa == ranking0_bv == ranking0_f64v == ranking0_org
@show gains0_pa == gains0_bv == gains0_f64v == gains0_org;

# %%
