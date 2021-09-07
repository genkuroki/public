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
function collatz(n)
    print(n)
    while n ∉ (1, 0, -1, -5, -17)
        n = iseven(n) ? n ÷ 2 : 3n + 1
        print("→", n)
    end
end

# %%
collatz(73)

# %%
collatz(-15)

# %%
collatz(-7)

# %%
collatz(-50)

# %%
collatz(big(2)^73 + 1)

# %%
@code_warntype collatz(73)

# %%
@code_warntype collatz(big(2)^73 + 1)

# %%
function collatz_seq(n, maxiters = 10^10)
    s = [n]
    l = 0
    while n > 1 && l < maxiters
        n = iseven(n) ? n ÷ 2 : 3n + 1
        push!(s, n)
        l += 1
    end
    s
end

collatz_seq(73)

# %%
function collatz_length(n, maxiters = 10^10)
    l = 0
    while n > 1 && l < maxiters
        n = iseven(n) ? n ÷ 2 : 3n + 1
        l += 1
    end
    l
end

collatz_length(73)

# %%
@code_warntype collatz_length(73, 10^10)

# %%
@code_native debuginfo=:none collatz_length(73, 10^10)

# %%
using Plots
n = 1:10^4
scatter(n, collatz_length.(n); label="", ms=1.5, msc=:auto, ma=0.5)

# %%
