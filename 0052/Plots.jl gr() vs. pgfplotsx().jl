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
#     display_name: Julia 1.11.3
#     language: julia
#     name: julia-1.11
# ---

# %%
using SpecialFunctions
using Plots
default(fmt=:png, legend=:topleft,
    legendfontsize=16, guidefontsize=16, tickfontsize=12, titlefontsize=18)
using LaTeXStrings

gr()
const mincho = "ipamp"
const gothic = "ipagp"

nanclamp(x, lo, hi) = lo ≤ x ≤ hi ? x : oftype(x, NaN)
stirling_approx(x) = x*log(x) - x - 0.5log(x) + 0.5log(2π)
xmin, xmax = 0, 10
ymin, ymax = -1, logabsgamma(xmax)[1]
plot(x -> nanclamp(logabsgamma(x)[1], ymin, ymax), xmin, xmax; label=L"y = \log\Gamma(x)")
plot!(x -> nanclamp(stirling_approx(x), ymin, ymax), xmin, xmax; label="Stirling近似", ls=:dash)
plot!(xtick=0:10, ytick=0:13)
plot!(xguide=L"x", yguide=L"y")
title!(L"Stirling近似: $\log\Gamma(x)\approx x\log x - x - \frac{1}{2}\log x + \frac{1}{2}\log 2\pi$")
annotate!(3, 3, text(raw"こんな感じによく近似されます。", :left, :red, 16, rotation=37, gothic))
plot!(fontfamily=mincho, size=(720, 480))


# %%
using SpecialFunctions
using Plots
default(fmt=:png, legend=:topleft,
    legendfontsize=16, guidefontsize=16, tickfontsize=12, titlefontsize=18)
using LaTeXStrings

gr()
const mincho = "ipamp"
const gothic = "ipagp"

nanclamp(x, lo, hi) = lo ≤ x ≤ hi ? x : oftype(x, NaN)
stirling_approx(x) = x*log(x) - x - 0.5log(x) + 0.5log(2π)
xmin, xmax = 0, 10
ymin, ymax = -1, logabsgamma(xmax)[1]
plot(x -> nanclamp(logabsgamma(x)[1], ymin, ymax), xmin, xmax; label=L"y = \log\Gamma(x)")
plot!(x -> nanclamp(stirling_approx(x), ymin, ymax), xmin, xmax; label="Stirling近似", ls=:dash)
plot!(xtick=0:10, ytick=0:13)
plot!(xguide=L"x", yguide=L"y")
title!(L"Stirling近似: $\log\Gamma(x)\approx x\log x - x - \frac{1}{2}\log x + \frac{1}{2}\log 2\pi$")
annotate!(3, 3, text(LaTeXString(raw"こんな感じによく近似されます。"), :left, :red, 16, rotation=37, gothic))
plot!(fontfamily=mincho, size=(720, 480))


# %%
using SpecialFunctions
using Plots
default(fmt=:png, legend=:topleft,
    legendfontsize=16, guidefontsize=16, tickfontsize=12, titlefontsize=18)
using LaTeXStrings

pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]

nanclamp(x, lo, hi) = lo ≤ x ≤ hi ? x : oftype(x, NaN)
stirling_approx(x) = x*log(x) - x - 0.5log(x) + 0.5log(2π)
xmin, xmax = 0, 10
ymin, ymax = -1, logabsgamma(xmax)[1]
plot(x -> nanclamp(logabsgamma(x)[1], ymin, ymax), xmin, xmax; label=L"y = \log\Gamma(x)")
plot!(x -> nanclamp(stirling_approx(x), ymin, ymax), xmin, xmax; label="Stirling近似", ls=:dash)
plot!(xtick=0:10, ytick=0:13)
plot!(xguide=L"x", yguide=raw"y")
title!(L"Stirling近似: $\log\Gamma(x)\approx x\log x - x - \frac{1}{2}\log x + \frac{1}{2}\log 2\pi$")
annotate!(3, 3, text(raw"\textbf{こんな感じによく近似されます。}", :left, :red, 16, rotation=37))


# %%
using SpecialFunctions
using Plots
default(fmt=:png, legend=:topleft,
    legendfontsize=16, guidefontsize=16, tickfontsize=12, titlefontsize=18)
using LaTeXStrings

pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]

nanclamp(x, lo, hi) = lo ≤ x ≤ hi ? x : oftype(x, NaN)
stirling_approx(x) = x*log(x) - x - 0.5log(x) + 0.5log(2π)
xmin, xmax = 0, 10
ymin, ymax = -1, logabsgamma(xmax)[1]
plot(x -> nanclamp(logabsgamma(x)[1], ymin, ymax), xmin, xmax; label=L"y = \log\Gamma(x)")
plot!(x -> nanclamp(stirling_approx(x), ymin, ymax), xmin, xmax; label=LaTeXString("Stirling近似"), ls=:dash)
plot!(xtick=0:10, ytick=0:13)
plot!(xguide=L"x", yguide=raw"y")
title!(L"Stirling近似: $\log\Gamma(x)\approx x\log x - x - \frac{1}{2}\log x + \frac{1}{2}\log 2\pi$")
annotate!(3, 3, text(LaTeXString(raw"\textbf{こんな感じによく近似されます。}"), :left, :red, 16, rotation=37))


# %%
using SpecialFunctions
using Plots
default(fmt=:png, legend=:topleft,
    legendfontsize=16, guidefontsize=16, tickfontsize=12, titlefontsize=18)
using LaTeXStrings

pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]
@eval Plots pgfx_sanitize_string(s::AbstractString) = s

nanclamp(x, lo, hi) = lo ≤ x ≤ hi ? x : oftype(x, NaN)
stirling_approx(x) = x*log(x) - x - 0.5log(x) + 0.5log(2π)
xmin, xmax = 0, 10
ymin, ymax = -1, logabsgamma(xmax)[1]
plot(x -> nanclamp(logabsgamma(x)[1], ymin, ymax), xmin, xmax; label=L"y = \log\Gamma(x)")
plot!(x -> nanclamp(stirling_approx(x), ymin, ymax), xmin, xmax; label="Stirling近似", ls=:dash)
plot!(xtick=0:10, ytick=0:13)
plot!(xguide=L"x", yguide=raw"y")
title!(L"Stirling近似: $\log\Gamma(x)\approx x\log x - x - \frac{1}{2}\log x + \frac{1}{2}\log 2\pi$")
annotate!(3, 3, text(raw"\textbf{こんな感じによく近似されます。}", :left, :red, 16, rotation=37))


# %%
