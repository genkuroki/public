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

# %% [markdown]
# https://discourse.julialang.org/t/why-isnt-a-simple-map-type-stable/67028
#
# https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured-1

# %%
function gen1(n)
    prevx = 0.0
    map(1:n) do _
        x = randn()
        r = prevx + x
        prevx = x
        r
    end
end

@code_warntype gen1(10)

# %%
function gen1_rev1(n)
    prevx = fill(0.0)
    map(1:n) do _
        x = randn()
        r = prevx[] + x
        prevx[] = x
        r
    end
end

@code_warntype gen1_rev1(10)

# %%
function gen1_rev2(n)
    prevx = Ref(0.0)
    map(1:n) do _
        x = randn()
        r = prevx[] + x
        prevx[] = x
        r
    end
end

@code_warntype gen1_rev2(10)

# %%
function gen1_rev3(n)
    prevx::Float64 = 0.0
    map(1:n) do _
        x = randn()
        r = prevx + x
        prevx = x
        r
    end
end

@code_warntype gen1_rev3(10)

# %%
