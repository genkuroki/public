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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# https://twitter.com/snap_tck/status/1586506933439434752

# %%
mrump(a, b) = 333.75b^6 + a^2*(11a^2*b^2 - b^6 - 121*b^4 - 2) + 5.5b^8
@show mrump(77617.0, 33096.0)
@show mrump(big"77617.0", big"33096.0")
@show mrump(big(77617), big(33096));

# %%
mrump_int(a, b) = (1335b^6 + 4a^2*(11a^2*b^2 - b^6 - 121*b^4 - 2) + 22b^8)/4
@show mrump_int(77617.0, 33096.0)
@show mrump_int(big"77617.0", big"33096.0")
@show mrump_int(big(77617), big(33096));

# %%
mrump_rat(a, b) = (1335//4)*b^6 + a^2*(11a^2*b^2 - b^6 - 121*b^4 - 2) + (22//4)*b^8
@show mrump_rat(77617.0, 33096.0)
@show mrump_rat(big"77617.0", big"33096.0")
@show mrump_rat(Int128(77617), Int128(33096))
@show mrump_rat(big(77617), big(33096));

# %%
precision(BigFloat)

# %%
setprecision(120) do
    mrump(big"77617.0", big"33096.0")
end

# %%
setprecision(124) do
    mrump(big"77617.0", big"33096.0")
end

# %%
@code_typed mrump(77617.0, 33096.0)

# %%
@code_typed mrump(big"77617.0", big"33096.0")

# %%
@code_typed mrump_rat(Int128(77617), Int128(33096))

# %%
