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

# %% [markdown]
# https://github.com/JuliaStats/StatsModels.jl/issues/220

# %% [markdown]
# ## Reproduce the slow-down (1)

# %%
using DataFrames
using GLM
using StatsBase: sample

n = 20
r = 10
x_symbols = [Symbol("x$i") for i in 1:n]
@time df = DataFrame(rand(100, n+1), [:y; x_symbols]);

# %%
x_vars = sample(x_symbols, r; replace=false)
@time F = term(:y) ~ sum(term(x) for x in x_vars)
@time cols = Tables.columntable(df)
@time mf = ModelFrame(F, cols, model=LinearModel)
@time mm = ModelMatrix(mf)
@time y = response(mf)
@time linmodel = fit(LinearModel, mm.m, y)
@time regmodel = StatsModels.TableRegressionModel(linmodel, mf, mm)

# %%
x_vars = sample(x_symbols, r; replace=false)
@time F = term(:y) ~ sum(term(x) for x in x_vars)
@time cols = Tables.columntable(df)
@time mf = ModelFrame(F, cols, model=LinearModel)
@time mm = ModelMatrix(mf)
@time y = response(mf)
@time linmodel = fit(LinearModel, mm.m, y)
@time regmodel = StatsModels.TableRegressionModel(linmodel, mf, mm)

# %% [markdown]
# ## Reproduce the slow-down (2)

# %%
using DataFrames
using GLM
using StatsBase: sample

n = 20
r = 10
x_symbols = [Symbol("x$i") for i in 1:n]
df = DataFrame(rand(100, n+1), [:y; x_symbols]);

result = []
for _ in 1:10
    x_vars = sample(x_symbols, r; replace=false)
    F = term(:y) ~ sum(term(x) for x in x_vars)
    @time regmodel = lm(F, df)
    push!(result, regmodel)
end
result

# %% [markdown]
# ## Solution

# %%
using DataFrames
using GLM
using StatsBase: sample

struct MyLinearModel{T, Y, X}
    linmodel::T
    y_var::Y
    x_vars::X
end

function my_lm(y_var, x_vars, df)
    y = df[!, y_var]
    X = [ones(nrow(df)) Matrix(df[!, x_vars])]
    linmodel = lm(X, y)
    MyLinearModel(linmodel, y_var, x_vars)
end

function Base.show(io::IO, mylm::MyLinearModel)
    linmodel, y_var, x_vars = getfield.(Ref(mylm), (:linmodel, :y_var, :x_vars))
    ct = coeftable(linmodel)
    ct.rownms .= string.((Symbol("(Intercept)"), x_vars...,))
    print(io, "Formula: ")
    print(io, y_var, " ~ 1")
    for x in x_vars print(io, " + ", x) end
    print(io, "\n\n")
    show(io, ct)
    print(io, "\n")
end

n = 20
r = 10
x_symbols = [Symbol("x$i") for i in 1:n]
df = DataFrame(rand(100, n+1), [:y; x_symbols]);

myresult = []
for _ in 1:10
    x_vars = sample(x_symbols, r; replace=false)
    @time mylinmodel = my_lm(:y, x_vars, df)
    push!(myresult, mylinmodel)
end
myresult

# %% [markdown]
# ## Analysis of the slow-down

# %%
x_vars = sample(x_symbols, r; replace=false)
@show x_vars
@time F = term(:y) ~ sum(term(x) for x in x_vars)
@time lm(F, df)

# %%
@time my_lm(:y, x_vars, df)

# %%
x_vars = sample(x_symbols, r; replace=false)
@time my_lm(:y, x_vars, df)

# %%
x_vars = sample(x_symbols, r; replace=false)
@time F = term(:y) ~ sum(term(x) for x in x_vars)
@time cols = Tables.columntable(df)
@time mf = ModelFrame(F, cols, model=LinearModel)
@time mm = ModelMatrix(mf)
@time y = response(mf)
@time linmodel = fit(LinearModel, mm.m, y)
@time regmodel = StatsModels.TableRegressionModel(linmodel, mf, mm);

# %%
typeof(cols)

# %%
typeof(mf)

# %%
typeof(mm) |> x -> (fieldnames(x), fieldtypes(x))

# %%
typeof(y)

# %%
typeof(linmodel)

# %%
typeof(regmodel)

# %%
regmodel.model == linmodel

# %%
regmodel.mf == mf

# %%
regmodel.mm == mm

# %%
@which lm(F, df)

# %%
@which fit(LinearModel, F, df)

# %%
@which fit(LinearModel, mm.m, y)

# %%
@code_warntype term(:y) ~ sum(term(x) for x in x_vars)

# %%
@code_warntype Tables.columntable(df)

# %%
@code_warntype ModelFrame(F, cols, model=LinearModel)

# %%
@code_warntype ModelMatrix(mf)

# %%
@code_warntype response(mf)

# %%
@code_warntype fit(LinearModel, mm.m, y)

# %%
@code_warntype StatsModels.TableRegressionModel(linmodel, mf, mm)

# %%
@code_warntype my_lm(:y, x_vars, df)

# %%
