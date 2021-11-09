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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# # Maxima.jl のインストールの仕方 (Windows編)
#
# 黒木玄
#
# 2018-01-03, 2019-06-08, 2019-12-27, 2021-11-09
#
# Old Version: https://nbviewer.org/gist/genkuroki/6d32a79a04770dc63d366c2f1a0d61a4
#
# (0) 以下にアクセスして色々読む.
#
# https://github.com/nsmith5/Maxima.jl
#
# https://nsmith5.github.io/Maxima.jl/latest/
#
# (1) Maxima のバイナリをインストールする.
#
# https://sourceforge.net/projects/maxima/files/Maxima-Windows/
#
# 注意 2019-12-27：これを書いている時点の最新版 maxima-clisp-sbcl-5.41.0a-win64.exe をインストールすると gnuplot で作画できなくなる. それより少し古い maxima-clisp-sbcl-5.38.1.exe であれば gnuplot で作画可能である.
#
# 注意 2021-11-09: maxima-5.45.1s-win64.exe をダウンロードしてこのノートブックを正常に実行できている。
#
# (2) Maximaがインストールされているディレクトリ以下を探して maxima.bat を見付け, maxima.bat がおいてあるディレクトリを環境変数 PATH に追加する.
#
# (3) julia で
#
# ```julia
# run('maxima')
# ```
#
# を実行して, julia から maxima.bat が起動できることを確認する.
#
# (4) julia で以下を実行する.
#
# ```julia
# using Pkg
# Pkg.add("Maxima")
# ```
#
# (5) 以下を実行できるかを確認する:
#
# ```julia
# using Maxima
# a = m"sin(%pi/3)"
# b = integrate(m"log(sin(x))", 'x')
# c = integrate(m"x*log(sin(x))", 'x')
# ```
#
# Windows が「アクセスを許可してよいか？」のような意味の警告を発して来た場合には「許可」するように答える.
#
# 追記 2021-11-09: 次の一行を追加しないと、 `display(x)` を手動で実行しないとエラーになってしまう。
#
# ```julia
# Base.convert(::Type{Any}, x::MExpr) = x
# ```

# %%
using Maxima
Base.convert(::Type{Any}, x::MExpr) = x

# %%
a = m"sin(%pi/3)"

# %%
b = integrate(m"log(sin(x))", 'x')

# %%
c = integrate(m"x*log(sin(x))", 'x')

# %%
d = factor(MExpr("a^2 - b^2"))

# %%
A = diff(m"log(f(x))", 'x', 4)

# %%
B = integrate(m"log(sin(x))", 'x', m"0", m"%pi/2")

# %%
C = integrate(m"1/sqrt(1+x+x^2)", 'x')

# %%
# gnuplot が立ち上がる
plot2d(m"sin(x)", x=(-4,4), title="Sine Wave")

# %%
# gnuplot が立ち上がる
plot3d(
m"[5*cos(x)*(cos(x/2)*cos(y)+sin(x/2)*sin(2*y)+3.0)-10.0,
-5*sin(x)*(cos(x/2)*cos(y)+sin(x/2)*sin(2*y)+3.0),
5*(-sin(x/2)*cos(y)+cos(x/2)*sin(2*y))], [x, -%pi, %pi], [y, -%pi, %pi]")

# %%
using Base64
showimg(mime, fn; tag="img") = open(fn) do f
    base64 = base64encode(f)
    display("text/html", """<$tag src="data:$mime;base64,$base64" />""")
end

showimg("image/png", "2021-11-09_00.png")
showimg("image/png", "2021-11-09_01.png")

# %%
