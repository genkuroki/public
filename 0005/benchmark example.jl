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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
function g!(x, y)
    for i in 1:length(x)
        x[i] *= 3
        y[i] += x[i]
    end
    x, y
end

function h!()
    for i in 1:length(x_global)
        x_global[i] *= 3
        y_global[i] += x_global[i]
    end
end

function f0()
    x = randn(100, 100)
    y = randn(100, 100)
    x, y
end

function f1()
    x = randn(100, 100)
    y = randn(100, 100)    
    g!(x, y)
end

function f2!(x = randn(100, 100), y = randn(100, 100))
    g!(x, y)
end

function f3!()
   g!(x_global, y_global)
end

function f4!()
   h!()
end

using BenchmarkTools

x_global = randn(100, 100)
y_global = randn(100, 100)

@show VERSION
@btime g!($x_global, $y_global)
@btime h!()
@btime f0()
@btime f1()
@btime f2!()
@btime f3!()
@btime f4!()

# %%
@code_warntype h!()

# %%
@code_warntype g!(x_global, y_global)

# %%
@code_typed f3!()

# %%
@code_typed g!(x_global, y_global)

# %%
@code_warntype f3!()

# %%
