# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %%
using Plots
function anim_stdratio(; n = 30)
    @gif for p̂ in [0.05:0.01:0.95; 0.95:-0.01:0.05]
        se = √(p̂*(1-p̂)/n)
        plot(p -> √(p*(1-p)/(p̂*(1-p̂))), max(0, p̂-3se), min(1, p̂+3se); label="")
        plot!(ylim=(0, 1.5))
        plot!(xguide="p", yguide="√(p(1-p)/(p̂(1-p̂)))")
        title!("p̂ = $p̂")
    end
end

# %%
anim_stdratio(n = 30)

# %%
anim_stdratio(n = 100)

# %%
anim_stdratio(n = 300)

# %%
anim_stdratio(n = 1000)

# %%
