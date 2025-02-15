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
# 文字列をLaTeXStringで十分にラップしたり、
# `@eval Plots pgfx_sanitize_string(s::AbstractString) = s` としていない場合には、
# {, } が \{, \} に変換されてしまい、「じ」「で」などの濁点の位置がおかしくなる。
# TeXのコードでは以下のようになってしまう。
# 
# title={\textbf\{$こ$$ん$$な$$感$$じ$$に$$日$$本$$語$$を$pgfplotsx()$で$$も$$使$$え$$ま$$す$!\}}

using Plots
default(fmt=:png)
using LaTeXStrings

pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]

plot(sin; label=L"正弦関数 $y = \sin\theta$", legendfontsize=15)
plot!(xguide=L"\theta", yguide=L"y", guidefontsize=14)
plot!(xtick=(-2π:π/2:2π, [L"-2\pi", L"-3\pi/2", L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),
    ytick=-1:0.2:1, tickfontsize=14)
plot!(ylim=(-1.2, 1.4), legend=:topleft)
annotate!(0, -0.3, text(L"正弦関数 $y = \sin\theta$", :left, 14, :red, rotation=45))
title!(raw"\textbf{こんな感じに日本語をpgfplotsx()でも使えます!}", titlefontsize=17)
plot!(tex_output_standalone=true)
savefig("pgfplotsx-ja-standalone-inappropriate.tex")
plot!()

# %%
; cat pgfplotsx-ja-standalone-inappropriate.tex

# %%
# 文字列をLaTeXStringで十分にラップすると日本語やLaTeXのコードを適切に使用できる.

using Plots
default(fmt=:png)
using LaTeXStrings

pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]

plot(sin; label=L"正弦関数 $y = \sin\theta$", legendfontsize=15)
plot!(xguide=L"\theta", yguide=L"y", guidefontsize=14)
plot!(xtick=(-2π:π/2:2π, [L"-2\pi", L"-3\pi/2", L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),
    ytick=-1:0.2:1, tickfontsize=14)
plot!(ylim=(-1.2, 1.4), legend=:topleft)
annotate!(0, -0.3, text(L"正弦関数 $y = \sin\theta$", :left, 14, :red, rotation=45))
title!(LaTeXString(raw"\textbf{こんな感じに日本語をpgfplotsx()でも使えます!}"), titlefontsize=17)
savefig("pgfplotsx-ja-fig.tex")
savefig("pgfplotsx-ja-fig.pdf")
plot!()

# %%
; cat pgfplotsx-ja-fig.tex

# %%
# そのままコンパイルできるLaTeXファイルの出力

using Plots
default(fmt=:png)
using LaTeXStrings

pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]

plot(sin; label=L"正弦関数 $y = \sin\theta$", legendfontsize=15)
plot!(xguide=L"\theta", yguide=L"y", guidefontsize=14)
plot!(xtick=(-2π:π/2:2π, [L"-2\pi", L"-3\pi/2", L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),
    ytick=-1:0.2:1, tickfontsize=14)
plot!(ylim=(-1.2, 1.4), legend=:topleft)
annotate!(0, -0.3, text(L"正弦関数 $y = \sin\theta$", :left, 14, :red, rotation=45))
title!(LaTeXString(raw"\textbf{こんな感じに日本語をpgfplotsx()でも使えます!}"), titlefontsize=17)
plot!(tex_output_standalone=true)
savefig("pgfplotsx-ja-standalone.tex")
plot!()

# %%
; cat pgfplotsx-ja-standalone.tex

# %%
# `@eval Plots pgfx_sanitize_string(s::AbstractString) = s` としてもよい。

using Plots
default(fmt=:png)
using LaTeXStrings

pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]
@eval Plots pgfx_sanitize_string(s::AbstractString) = s

# @eval Plots function original_pgfx_sanitize_string(s::AbstractString)
#     # regular latex text with the following special characters won't compile if not sanitized (escaped)
#     sanitized = replace(s, r"\\?([#%_&\{\}\$])" => s"\\\1")
#     map(collect(sanitized)) do c
#         if isascii(c)
#             c
#         else
#             Latexify.latexify(c; parse = false)
#         end
#     end |> join |> LaTeXString
# end
# @eval Plots pgfx_sanitize_string(s::AbstractString) = original_pgfx_sanitize_string(s)

plot(sin; label=L"正弦関数 $y = \sin\theta$", legendfontsize=15)
plot!(xguide=L"\theta", yguide=L"y", guidefontsize=14)
plot!(xtick=(-2π:π/2:2π, [L"-2\pi", L"-3\pi/2", L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),
    ytick=-1:0.2:1, tickfontsize=14)
plot!(ylim=(-1.2, 1.4), legend=:topleft)
annotate!(0, -0.3, text(L"正弦関数 $y = \sin\theta$", :left, 14, :red, rotation=45))
title!(raw"\textbf{こんな感じに日本語をpgfplotsx()でも使えます!}", titlefontsize=17)
plot!()

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
