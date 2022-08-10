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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
function f(a, b)
    s = zero(a+b)
    for i in a:b
        s += i
    end
    s
end

# %%
function f_annotated(a::Int, b::Int)
    s = 0
    for i in a:b
        s += i
    end
    s
end

# %%
f(10, 20)

# %%
f_annotated(10, 20)

# %%
@code_native debuginfo=:none f(10, 20)

# %%
@code_native debuginfo=:none f_annotated(10, 20)

# %% tags=[]
f(10.0, 20.0)

# %% tags=[]
f_annotated(10.0, 20.0)

# %%
F(a, b) = zero(a + b)

# %%
F(10, 20)

# %%
F(10.0, 20)

# %%
@code_llvm F(10, 20)

# %%
@code_llvm F(10.0, 20)

# %%
G(a, b) = zero(eltype(a:b))

# %%
@code_llvm F(10.0, 20)

# %%
function f1(a, b)
    s = zero(eltype(a:b))
    for i in a:b
        s += i
    end
    s
end

# %%
@code_native debuginfo=:none f1(10, 20)

# %%
using Dates

# %%
x = Day(7)

# %%
zero(x)

# %%
one(x)

# %%
oneunit(x)

# %%
