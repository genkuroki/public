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
# Straightforward translation of https://github.com/rmeinl/apricot-julia/blob/5f130f846f8b7f93bb4429e2b182f0765a61035c/notebooks/python_reimpl.ipynb

# %%
using Seaborn
using ScikitLearn: @sk_import
@sk_import datasets: fetch_covtype
using Random
using StatsBase

# %%
digits_data = fetch_covtype()

# %%
X_digits = permutedims(abs.(digits_data["data"]))
summary(X_digits)

# %%
"""`calculate_gains!(X, gains, current_values, idxs, current_concave_values_sum)` mutate `gains` only"""
function calculate_gains!(X, gains, current_values, idxs, current_concave_values_sum)
    Threads.@threads for i in eachindex(idxs)
        @inbounds idx = idxs[i]
        @inbounds gains[i] = sum(j -> sqrt(current_values[j] + X[j, idx]), axes(X, 1))
    end
    gains .-= current_concave_values_sum
end

@doc calculate_gains!

# %%
function fit(X, k; calculate_gains! = calculate_gains!)
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

        popat!(idxs, best_idx)
    end
    return ranking, total_gains
end

# %%
k = 1000

# %%
@time ranking0, gains0 = fit(X_digits, k; calculate_gains! = calculate_gains!);

# %%
@time ranking0, gains0 = fit(X_digits, k; calculate_gains! = calculate_gains!);

# %%
tic = time()
ranking0, gains0 = fit(X_digits, k; calculate_gains! = calculate_gains!)
toc0 = time() - tic
toc0

# %%
@time begin
idxs = sample(axes(X_digits, 2), k; replace=false)
X_subset = X_digits[:, idxs]
gains1 = cumsum(X_subset; dims=2)
gains1 = vec(sum(sqrt, gains1; dims=1))
end;

# %%
@time begin
idxs = sample(axes(X_digits, 2), k; replace=false)
X_subset = X_digits[:, idxs]
gains1 = cumsum(X_subset; dims=2)
gains1 = vec(sum(sqrt, gains1; dims=1))
end;

# %%
tic = time()
idxs = sample(axes(X_digits, 2), k; replace=false)
X_subset = X_digits[:, idxs]
gains1 = cumsum(X_subset; dims=2)
gains1 = vec(sum(sqrt, gains1; dims=1))
toc1 = time() - tic
toc1

# %%
plt.figure(figsize=(9, 4.5))

plt.subplot(121)
plt.plot(cumsum(gains0), label="Naive")
plt.plot(gains1, label="Random")
plt.ylabel("F(S)")
plt.xlabel("Subset Size")
plt.legend()
plt.grid(lw=0.3)

plt.subplot(122)
plt.bar(1:2, [toc0,  toc1])
plt.ylabel("Time (s)")
plt.xticks(1:2, ["Naive", "Random"], rotation=90)
plt.grid(lw=0.3)
plt.title("Julia 1.6.2")

plt.tight_layout()
#plt.show()

# %%
