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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %%
; ls -l hand-written-greek-letters.jpg

# %%
; magick hand-written-greek-letters.jpg hand-written-greek-letters.ppm

# %%
; ls -l hand-written-greek-letters.ppm

# %%
; potrace -s hand-written-greek-letters.ppm -o hand-written-greek-letters.svg

# %%
; ls -l hand-written-greek-letters.svg

# %%
; potrace -b pdf hand-written-greek-letters.ppm -o hand-written-greek-letters.pdf

# %%
; ls -l hand-written-greek-letters.pdf

# %%
