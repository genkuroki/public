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
function f()
    local i
    for outer i in 1:10^6
        i > 5*10^5 && break
    end
    i
end

# %%
function g()
    i = undef
    for outer i in 1:10^6
        i > 5*10^5 && break
    end
    i
end

# %%
f()

# %%
@code_warntype f()

# %%
@code_warntype g()

# %%
using BenchmarkTools

# %%
@btime f()

# %%
@btime g()

# %%
@code_native debuginfo=:none f()

# %%
@code_native debuginfo=:none g()

# %%
