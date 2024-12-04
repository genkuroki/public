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
#     display_name: Julia 1.11.1
#     language: julia
#     name: julia-1.11
# ---

# %%
n =

println(3)

# %%
println(n)

# %%
:(
n =

println(3)
)

# %%
:(
n =

println(3)
) |> eval

# %%
:(
n =

println(3)
) |> dump

# %%
:(
n =

println(3)
) |> Meta.show_sexpr

# %%
x =
1 +
2 *
3

# %%
:(
x =
1 +
2 *
3
)

# %%
:(
x =
1 +
2 *
3
) |> eval

# %%
:(
x =
1 +
2 *
3
) |>dump

# %%
:(
x =
1 +
2 *
3
) |> Meta.show_sexpr

# %%
where where where where where where where

# %%
:(where where where where where where where)

# %%
:(where where where where where where where) |> dump

# %%
:(where where where where where where where) |> Meta.show_sexpr

# %%
