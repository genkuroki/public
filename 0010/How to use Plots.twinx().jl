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
using Plots.PlotMeasures
pyplot(fmt = :svg)
bar(rand(10), c=:white, right_margin=1cm)
plot!(twinx(), rand(10), label="v1")
plot!(twinx(), rand(10) * 10, label="v2")

# %%
using Plots
using Plots.PlotMeasures
pyplot(fmt = :svg)
bar(rand(10), c=:white, right_margin=1cm)
plot!(twinx(), [rand(10) rand(10) * 10], label=["v1" "v2"])

# %%
using Plots
using Plots.PlotMeasures
pyplot(fmt = :svg)
plot(rand(10), label="v1", right_margin=1cm)
plot!(rand(10) * 10, label="v2")
bar!(twinx(), rand(10), label="", fillalpha=0)

# %%
using Plots
using Plots.PlotMeasures
pyplot(fmt = :svg)
bar(rand(10), c=:white, right_margin=1cm)
ax2 = twinx()
plot!(ax2, rand(10), label="v1")
plot!(ax2, rand(10) * 10, label="v2")

# %%
