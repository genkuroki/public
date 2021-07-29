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
module A

function f(x)
    a = x^2
    eval(:(a + 1))
end

end

A.f(3)

# %%
module A1

function f(x)
    global a = x^2
    eval(:(a + 1))
end

end

A1.f(3)

# %%
module B

function f(x)
    a = x^2
    eval(:($a + 1))
end

end

B.f(3)

# %%
@code_warntype B.f(3)

# %%
module C

function f(x)
    a = x^2
    eval(Expr(:call, :+, a, 1))
end

end

C.f(3)

# %%
@code_warntype C.f(3)

# %%
module D

function f(x)
    a = x^2
    a + 1
end

end

D.f(3)

# %%
@code_warntype D.f(3)

# %%
z = 9
:($z + 1) == Expr(:call, :+, z, 1)

# %%
