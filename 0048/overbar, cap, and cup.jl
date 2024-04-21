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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# https://x.com/dannchu/status/1780507090496409819

# %%
U = [1,2,3,4,5,6,7,8,9]
A = [1,2,3,5,7]
B = [1,3,4,5]
A̅ = setdiff(U,A) # A̅ is A\overbar + tab-key
B̅ = setdiff(U,B) # B̅ is B\overbar + tab-key

println("(1) A̅ =" , A̅ )
println("(2) B̅ =" , B̅)
println("(3) A̅ ∩ B̅ =", A̅ ∩ B̅ ) # ∩ is \cap + tab-key
println("(4) A̅ ∪ B̅ =", A̅ ∪ B̅ ) # ∪ is \cup + tab-key
println("(5) A̅ ∩ B =", A̅ ∩ B ) 
println("(6) A ∩ B̅ =", A ∩ B̅ ) 
println("(7) A ∪ B̅ =", A ∪ B̅ )  
println("(8) A̅ ∪ B =", A̅ ∪ B )

# %%
?A̅

# %%
?∩

# %%
?∪

# %%
