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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using Plots
function res_plot(n; label="n = $n", kwargs...)
    plot(rand(n); label, kwargs...)
end
res_plot(10)

# %%
res_plot(10; label="rand(10)", seriestype=:bar)

# %%
function plot_MxN(M, N;
        size = (720, 360),
        layout = (M, N),
        labelfunc = (i, j) -> "($i, $j)",
        kwargs...
    )
    PP = []
    for i in 1:M, j in 1:N
        P = plot(rand(10); label = labelfunc(i, j))
        push!(PP, P)
    end
    plot(PP...; size, layout, kwargs...)
end

plot_MxN(3, 4; tickfontsize = 5)

# %%
plot_MxN(3, 4; tickfontsize = 5, labelfunc = (i, j) -> "\$P_{$i,$j}\$")

# %%
