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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
a = 1:9
a*a'

# %%
a = 2:15
a*a'

# %%
a.∉(a*a',)

# %%
a[a.∉(a*a',)]

# %%
2:100|>a->a[a.∉(a*a',)]

# %%
:(2:100|>a->a[a.∉(a*a',)])

# %%
quote
    2:100|>a->a[a.∉(a*a',)]
end

# %%
quote
    y = x < 0 ? -x : x
end

# %%
quote
    y = if x < 0; -x else x end
end

# %%
@code_warntype (a->a[a.∉(a*a',)])(2:100)

# %%
@code_typed (a->a[a.∉(a*a',)])(2:100)

# %%
