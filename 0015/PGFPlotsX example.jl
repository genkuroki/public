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
;ls D:\\texlive\\2021\\bin\\win32\\lualatex.exe

# %%
using PGFPlotsX

p1 = range(0, 360; length = 100)' # parameter 1
p2 = range(0, 360; length = 70)   # parameter 2
w1 = @. sind(3*(p1 + 2 * p2)) + 1.25 # wave 1
w2 = @. 6 + w1 * cosd(p1)            # wave 2
x = vec(@. w2 * cosd(p2))
y = vec(@. w2 * sind(p2))
z = vec(@. w1 * sind(p1))
logo_colors = [(77,100,174), (57,151,79), (255,255,255), (146,89,163), (202,60,50)]

@pgf Axis(
    {
        axis_equal,
        axis_lines = "none"
    },
    Plot3(
        {
            surf,
            z_buffer = "sort",
            colormap = "{Julia}{$(join(map(c -> "rgb255=$c", logo_colors), ", "))}",
            "mesh/rows" = length(p1)
        },
        Table(x, y, z)))

# %%
