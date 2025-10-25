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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
function f(n=10^8)
    s = 0.0
    for k in 1:n
        s += sin(k)/k
    end
    1 + 2s
end

@time f()
@time f()
@time f()

# %%
using LoopVectorization

function g(n=10^8)
    s = 0.0
    @turbo for k in 1:n
        s += sin(k)/k
    end
    1 + 2s
end

@time g()
@time g()
@time g()

# %%
??@turbo

# %%
