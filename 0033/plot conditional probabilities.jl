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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10)
safediv(x, y) = x==0 ? x : y==Inf ? zero(y) : x/y

# %% [markdown]
# https://twitter.com/03t_ms/status/1545720479868022784
#
# <img src="https://github.com/genkuroki/public/raw/main/0033/FXOAk2cUcAArHXd.jpg" width=70%>

# %%
RWB = [
    10  5  3
     8  4  6
     4  3  5
]
Z = [1//3, 1//3, 1//3] .* RWB .// sum(RWB; dims=2)

# %%
sum(Z; dims=2)

# %%
sum(Z)

# %%
display(Z)

X = reverse(Z; dims=1)
xtick = (1:3, ["R", "W", "B"])
ytick = (1:3, reverse(["A", "B", "C"]))
ann = [(j, i, text(X[i,j], 16, :cyan, :center))
            for i in 1:3 for j in 1:3]
O1 = heatmap(X; xtick, ytick)
annotate!(ann)
title!("joint probability distribution")

# %%
sum(Z; dims=1)

# %%
W = Z .// sum(Z; dims=1)

# %%
sum(W; dims=1)

# %%
display(W)

X = reverse(W; dims=1)
xtick = (1:3, ["red", "white", "blue"])
ytick = (1:3, reverse(["A", "B", "C"]))
ann = [(j, i, text(X[i,j], 16, :cyan, :center))
            for i in 1:3 for j in 1:3]
O2 = heatmap(X; xtick, ytick)
annotate!(ann)
title!("conditional probability distribution")

# %%
display(sum(Z; dims=1))
plot(O1, O2; size=(400, 500), layout=(2, 1))

# %%
for k in 1:3
    @eval @show $(Z[k,:]) / $(sum(Z; dims=1)[k])
end

# %%
n = 10^6
x = randn(n)
e = randn(n)
y = @. 1 + x + e

m = 10^4
xbin = -3:0.25:3
ybin = -2:0.25:4

P0 = scatter(x[1:m], y[1:m]; label="y = 1 + x + e", ms=2, ma=0.5, msw=0)
plot!(x -> 1 + x; label="y = 1 + x")
plot!(legend=:topleft, colorbar=false)
plot!(xlim=extrema(xbin), ylim=extrema(ybin))
title!("data")

# %%
h = Plots._make_hist((x, y), (xbin, ybin); normed=true)

P1 = plot(h)
plot!(x -> 1 + x; label="y = 1 + x", c=:cyan)
plot!(legend=:topleft, colorbar=false)
plot!(xlim=extrema(xbin), ylim=extrema(ybin))
title!("joint probability distribution")

# %%
H = deepcopy(h)
H.weights = safediv.(h.weights, sum(h.weights; dims=2))

P2 = plot(H)
plot!(x -> 1 + x; label="y = 1 + x", c=:cyan)
plot!(legend=:topleft, colorbar=false)
plot!(xlim=extrema(xbin), ylim=extrema(ybin))
title!("conditional probability distribution")

# %%
plot(P0, P1, P0, P2; size=(800, 600), layout=(2, 2))

# %%
