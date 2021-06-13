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
#     display_name: Julia 1.7.0-DEV
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# https://discourse.julialang.org/t/proper-modules-in-the-future/62646/112

# %%
dispcode(c, l = "julia") = display("text/markdown", "```$l\n$c\n```")

# %%
MyCalc_jl = raw"""
calculate(x) = x * subcalc(x)
"""
dispcode(MyCalc_jl)

# %%
subcalc(x) = x^2
include_string(Main, MyCalc_jl)
calculate(2)

# %%
MyCalc2_jl = raw"""
module MyCalc2
calculate2(x) = x * subcalc2(x)
end
"""
dispcode(MyCalc2_jl)

# %%
subcalc2(x) = x^2
include_string(Main, MyCalc2_jl)
using .MyCalc2: calculate2
calculate2(2)

# %% tags=[]
MyCalc3_jl = raw"""
module MyCalc3
calculate3(x)= x * subcalc3(x)
subcalc3(x) = x^3
end
"""
dispcode(MyCalc3_jl)

# %%
include_string(Main, MyCalc3_jl)
using .MyCalc3: calculate3
calculate3(2)

# %%
MyCalc4_jl = raw"""
module MyCalc4
calculate4(x) = x * subcalc4(x)
using Main: subcalc4
end
"""
dispcode(MyCalc4_jl)

# %%
subcalc4(x) = x^4
include_string(Main, MyCalc4_jl)
using .MyCalc4: calculate4
calculate4(2)

# %%
MyCalc5_jl = raw"""
module MyCalc5
calculate5(subcalc5, x) = x * subcalc5(x)
end
"""
dispcode(MyCalc5_jl)

# %%
include_string(Main, MyCalc5_jl)
using .MyCalc5: calculate5
calculate5(x -> x^5, 2)

# %%
