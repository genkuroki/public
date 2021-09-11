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

# %%
(a->a[a.∉(a*a',)])(2:100)

# %%
function supertypes(x)
    T = typeof(x)
    t = DataType[]
    while T != Any
        push!(t, T)
        T = supertype(T)
    end
    t
end

# %%
a = 2:100

# %%
supertypes(a)

# %%
a'

# %%
supertypes(a')

# %%
a*a'

# %%
supertypes(a*a')

# %%
a .∉ (a*a',)

# %%
supertypes(a .∉ (a*a',))

# %%
a[a .∉ (a*a',)]

# %%
supertypes(a[a .∉ (a*a',)])

# %%
