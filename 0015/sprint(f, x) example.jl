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
function show_Hello(io::IO, x)
    print(io, "Hello, ")
    show(io, x)
    print(io, '!')
end
show_Hello(x) = show_Hello(stdout, x)

show_Hello("Julia")

# %%
str = sprint(show_Hello, "Julia")

# %%
print(str)

# %%
?sprint

# %%
methods(sprint)

# %%
