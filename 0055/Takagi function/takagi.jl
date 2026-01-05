# -*- coding: utf-8 -*-
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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
using Plots
default(fmt=:png, legend_font_halign=:left)
pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]
@eval Plots pgfx_sanitize_string(s::AbstractString) = s

function takagi(x; a=1/2, b=2.0, niters=20)
    y = zero(float(x))
    for n in 0:niters
        y += a^n * abs(b^n * x - round(b^n * x))
    end
    y
end

plot(takagi, 0, 1; label="")
plot!(; xtick=0:0.1:1, ytuck=0:0.1:1)
plot!(; xguide=raw"$x$", guidefontsize=16)
plot!(; yguide=raw"$\displaystyle
    T(x)=\sum_{n=0}^\infty\frac{s(2^n x)}{2^n}$")
title!("高木関数"; titlefontsize=20)
plot!(; tex_output_standalone=true)
savefig("takagi_standalone.tex")
plot!()

# %%
; lualatex takagi_standalone.tex

# %%
