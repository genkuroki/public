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
#     display_name: Julia 1.7.0-DEV
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# See https://twitter.com/taketo1024/status/1400990265364549634
#
# A straightforward translation in Swift to Julia of
# * https://gist.github.com/taketo1024/e1e1b92b0c5c9c441bd870476356204e

# %%
step = 10

function forAnySmall(p)
    for i in 0:step-1
        ε = (step - i) / step
        if !p(ε)
            return false
        end
    end
    return true
end

function existsSmall(p)
    for i in 0:step-1
        δ = (step - i) / step
        if p(δ) 
            return true
        end
    end
    return false
end

function forAnyInRange((from, to), p)
    for _ in 0:step-1
        x = from + rand() * (to - from)
        if !p(x)
            return false
        end
    end
    return true
end

function isContinuous(f, at, a)
    forAnySmall(ε ->
        existsSmall(δ ->
            forAnyInRange((a - δ, a + δ), x ->
                abs(f(x) - f(a)) < ε
            )
        )
    )
end

println(isContinuous(x -> x^2, :at, 0))
# true: x^2 is continuous at x = 0
println(isContinuous(floor, :at, 1))
# false: floor func is not continuous at x = 0

# %%
