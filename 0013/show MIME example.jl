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
A = rand(0:9, 3, 4)

# %%
display(A)

# %%
print("A = "); display(A)

# %%
print("A = "); flush(stdout); display(A)

# %%
print("A = "); show(stdout, MIME"text/plain"(), A)

# %%
print("A = "); show(stdout, MIME("text/plain"), A)

# %%
?MIME

# %%
