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
#     display_name: Julia 1.10.1
#     language: julia
#     name: julia-1.10
# ---

# %%
using CodeTracking
using Markdown
function mdcode(s; lang="julia")
    """
    ```$lang
    $s
    ```
    """ |> Markdown.parse
end

# %%
mdcode(@code_string sum(1:4))

# %%
A = rand(2,2)
mdcode(@code_string A * A)

# %%
