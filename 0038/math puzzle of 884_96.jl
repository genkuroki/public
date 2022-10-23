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

# %% [markdown]
# * https://twitter.com/884_96/status/1582718783344386054
# * https://twitter.com/dannchu/status/1584026377787146241

# %%
using Combinatorics
["$ア$イ × $ウ$エ$オ = $カ$キ$ク$ケ"
    for (ア,イ,ウ,エ,オ,カ,キ,ク,ケ) in permutations(1:9)
        if (10ア+イ) * (100ウ+10エ+オ) == 1000カ+100キ+10ク+ケ]

# %%
