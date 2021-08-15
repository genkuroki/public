# -*- coding: utf-8 -*-
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
# https://discourse.julialang.org/t/julia-beginner-from-python-numba-outperforms-julia-in-rewrite-any-tips-to-improve-performance/66414

# %%
using BenchmarkTools
#using Distributed
#using Plots
#using Profile
#using ProfileVega
#using PyCall
using ScikitLearn
#using Traceur
using TimerOutputs

# %%
@sk_import datasets: (fetch_covtype)

# %%
digits_data = fetch_covtype();

# %%
X_digits = abs.(digits_data["data"]);
X_digits = transpose(X_digits);

# %%
d, n = size(X_digits)

# %%
struct NaiveGreedy
    X::Matrix{Float64}
end

struct Result
    ranking::Vector{Int32}
    gains::Vector{Float64}
end

# %%
function lexsort(a, b, rev=false) 
    idxs = sortperm(a, alg=MergeSort, rev=rev)
    return idxs[sortperm(b[idxs], alg=MergeSort, rev=rev)]
end

function get_gains!(X, current_values, idxs, gains)
    @inbounds Threads.@threads for i in eachindex(idxs)
        s = 0.0
        for j in eachindex(current_values)
            s += @fastmath sqrt(current_values[j] + X[j, idxs[i]])
        end
        gains[i] = s
    end
end

function calculate_gains!(X, gains, current_values, idxs, current_concave_value_sum)
    get_gains!(X, current_values, idxs, gains)
    
    gains .-= current_concave_value_sum
    return gains
end

# %%
function fit(optimizer::NaiveGreedy, k, sample_cost)
    @timeit to "intro" begin
        d, n = size(optimizer.X)
        
        cost = 0.0

        ranking = Int32[]
        total_gains = Float64[]

        mask = zeros(Int8, n)
        current_values = zeros(Float64, d)
        current_concave_values = sqrt.(current_values)
        current_concave_values_sum = sum(current_concave_values)

        idxs = 1:n
    end

    gains = zeros(Float64, size(idxs)[1])
    @timeit to "while loop" begin
        while cost < k
            gains = @timeit to "calc_gains" calculate_gains!(optimizer.X, gains, current_values, idxs, current_concave_values_sum)
            
            if sample_cost != nothing
                gains ./= sample_cost[idxs]
                idx_idxs = @timeit to "lexsort" lexsort(gains, 1:size(gains)[1], true)

                @timeit to "select_idx" begin
                    for i in 1:size(idx_idxs)[1]
                        global idx = idx_idxs[i]
                        global best_idx = idxs[idx]
                        if cost + sample_cost[best_idx] <= k
                            break
                        end
                    end
                end
                curr_cost = sample_cost[best_idx]
            else
                global idx = argmax(gains)
                global best_idx = idxs[idx]
                curr_cost = 1.
            end
            
            if cost + curr_cost > k
                break
            end
                
            @timeit to "select_idx" begin
                cost += curr_cost
                # Calculate gains
                gain = gains[idx] * curr_cost
            end
            
            @timeit to "select_next" begin
                # Select next
                current_values += view(optimizer.X, :, best_idx)
                current_concave_values .= sqrt.(current_values)
                current_concave_values_sum = sum(current_concave_values)

                push!(ranking, best_idx)
                push!(total_gains, gain)

                mask[best_idx] = 1
                idxs = findall(==(0), mask)
            end
        end
    end
    return Result(ranking, total_gains)
end

# %%
const to = TimerOutput()
k = 1000
opt1 = NaiveGreedy(X_digits)
res1 = fit(opt1, k, nothing)
to

# %%
function fit(optimizer::NaiveGreedy, k, sample_cost)
    @timeit to "intro" begin
        d, n = size(optimizer.X)
        
        cost = 0.0

        ranking = Int32[]
        total_gains = Float64[]

        #mask = zeros(Int8, n)
        current_values = zeros(Float64, d)
        current_concave_values = sqrt.(current_values)
        current_concave_values_sum = sum(current_concave_values)

        #idxs = 1:n
        idxs = collect(1:n)
    end

    gains = zeros(Float64, size(idxs)[1])
    @timeit to "while loop" begin
        while cost < k
            gains = @timeit to "calc_gains" calculate_gains!(optimizer.X, gains, current_values, idxs, current_concave_values_sum)
            
            if sample_cost != nothing
                gains ./= sample_cost[idxs]
                idx_idxs = @timeit to "lexsort" lexsort(gains, 1:size(gains)[1], true)

                @timeit to "select_idx" begin
                    for i in 1:size(idx_idxs)[1]
                        global idx = idx_idxs[i]
                        global best_idx = idxs[idx]
                        if cost + sample_cost[best_idx] <= k
                            break
                        end
                    end
                end
                curr_cost = sample_cost[best_idx]
            else
                global idx = argmax(gains)
                global best_idx = idxs[idx]
                curr_cost = 1.
            end
            
            if cost + curr_cost > k
                break
            end
                
            @timeit to "select_idx" begin
                cost += curr_cost
                # Calculate gains
                gain = gains[idx] * curr_cost
            end
            
            @timeit to "select_next" begin
                # Select next
                #current_values += view(optimizer.X, :, best_idx)
                #current_concave_values .= sqrt.(current_values)
                #current_concave_values_sum = sum(current_concave_values)
                current_values .+= view(optimizer.X, :, best_idx)
                current_concave_values_sum = sum(sqrt, current_values)

                push!(ranking, best_idx)
                push!(total_gains, gain)

                #mask[best_idx] = 1
                #idxs = findall(==(0), mask)
                popat!(idxs, findfirst(==(best_idx), idxs))
            end
        end
    end
    return Result(ranking, total_gains)
end

# %%
const to = TimerOutput()
k = 1000
opt1rev = NaiveGreedy(X_digits)
res1rev = fit(opt1rev, k, nothing)
to

# %%
@show res1rev.ranking == res1.ranking
@show res1rev.gains ≈ res1.gains;

# %%
