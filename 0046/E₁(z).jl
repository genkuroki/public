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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # E₁(z)の教訓
#
# * 黒木玄
# * 2020-09-28, 2024-01-02 (rerun with Julia v1.10.0)
# $
# \newcommand\real{\operatorname{Re}}
# \newcommand\imag{\operatorname{Im}}
# $

# %% [markdown] toc=true
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#文脈" data-toc-modified-id="文脈-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>文脈</a></span><ul class="toc-item"><li><span><a href="#MITでの講義の宿題の模範解答" data-toc-modified-id="MITでの講義の宿題の模範解答-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>MITでの講義の宿題の模範解答</a></span></li><li><span><a href="#これは何を意味しているのか？" data-toc-modified-id="これは何を意味しているのか？-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>これは何を意味しているのか？</a></span></li><li><span><a href="#最適化の中身の紹介" data-toc-modified-id="最適化の中身の紹介-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>最適化の中身の紹介</a></span></li><li><span><a href="#コードの自動生成機能の重要性" data-toc-modified-id="コードの自動生成機能の重要性-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>コードの自動生成機能の重要性</a></span></li><li><span><a href="#教訓" data-toc-modified-id="教訓-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>教訓</a></span></li><li><span><a href="#「怠惰な人間」という巨大因子を無視してはいけない" data-toc-modified-id="「怠惰な人間」という巨大因子を無視してはいけない-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span><strong>「怠惰な人間」という巨大因子</strong>を無視してはいけない</a></span></li><li><span><a href="#uncorrelated氏による連分数のベンチマークテスト" data-toc-modified-id="uncorrelated氏による連分数のベンチマークテスト-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>uncorrelated氏による連分数のベンチマークテスト</a></span></li><li><span><a href="#補足と感謝" data-toc-modified-id="補足と感謝-1.8"><span class="toc-item-num">1.8&nbsp;&nbsp;</span>補足と感謝</a></span></li></ul></li><li><span><a href="#MITでの講義の宿題の解答のコード" data-toc-modified-id="MITでの講義の宿題の解答のコード-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>MITでの講義の宿題の解答のコード</a></span><ul class="toc-item"><li><span><a href="#連分数で計算する部分の解説" data-toc-modified-id="連分数で計算する部分の解説-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>連分数で計算する部分の解説</a></span></li></ul></li><li><span><a href="#uncorrelated氏のコード" data-toc-modified-id="uncorrelated氏のコード-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>uncorrelated氏のコード</a></span><ul class="toc-item"><li><span><a href="#最初のバージョンでは-reltol=1e-16-になっていた！(笑)" data-toc-modified-id="最初のバージョンでは-reltol=1e-16-になっていた！(笑)-3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>最初のバージョンでは reltol=1e-16 になっていた！(笑)</a></span></li><li><span><a href="#連分数展開が効率的でない領域で連分数のコードをテストしている！" data-toc-modified-id="連分数展開が効率的でない領域で連分数のコードをテストしている！-3.2"><span class="toc-item-num">3.2&nbsp;&nbsp;</span>連分数展開が効率的でない領域で連分数のコードをテストしている！</a></span></li><li><span><a href="#実践的には決して使われることがない遅いアルゴリズムで比較している！" data-toc-modified-id="実践的には決して使われることがない遅いアルゴリズムで比較している！-3.3"><span class="toc-item-num">3.3&nbsp;&nbsp;</span>実践的には決して使われることがない遅いアルゴリズムで比較している！</a></span></li></ul></li><li><span><a href="#uncorrelated氏の計算法がどれだけ遅くなっているか" data-toc-modified-id="uncorrelated氏の計算法がどれだけ遅くなっているか-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>uncorrelated氏の計算法がどれだけ遅くなっているか</a></span><ul class="toc-item"><li><span><a href="#素朴な比較" data-toc-modified-id="素朴な比較-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>素朴な比較</a></span></li><li><span><a href="#連分数展開が有効な領域での比較" data-toc-modified-id="連分数展開が有効な領域での比較-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>連分数展開が有効な領域での比較</a></span></li><li><span><a href="#結論:-論外！" data-toc-modified-id="結論:-論外！-4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>結論: 論外！</a></span></li></ul></li></ul></div>

# %%
using SymPy
using BenchmarkTools
using Plots
default(fmt=:png)
versioninfo()

# %%
# Override
# https://github.com/jverzani/SymPyCore.jl/blob/main/src/SymPy/show_sympy.jl#L31-L34
@eval SymPy begin
function Base.show(io::IO,  ::MIME"text/latex", x::SymbolicObject)
    out = _sympy_.latex(↓(x), mode="inline",fold_short_frac=false)
    out = replace(out, r"\\frac{"=>"\\dfrac{")
    print(io, string(out))
end
end

# %% [markdown]
# ## 文脈

# %% [markdown]
# ### MITでの講義の宿題の模範解答
#
# まずはMITでの講義の宿題の模範解答
#
# * https://nbviewer.jupyter.org/github/stevengj/18S096/blob/iap2017/pset3/pset3-solutions.ipynb
#
# を見て欲しい.  そこでは, 指数積分函数 E₁(z) のJuliaによる実装によって, scipyで採用されているFortranで書かれてたライブラリの5～6倍の速さを実現している.

# %% [markdown]
# ### これは何を意味しているのか？
#
# これを「Juliaのコンパイラの方がFortranよりも優れた最適化を行うことを意味する」と**誤解**してはいけない！
#
# ある程度以上の性能を持つコンパイラの最適化のレベルはオーダー的にはそう変わらない. 2倍の違いはないと思ってよい.
#
# それではどうしてJuliaで実装された指数積分函数 E₁(z) の方が Fortran で書かれたライブラリの5～6倍の速さで計算できてしまうのだろうか?
#
# Julia側が5～6倍も速くなっているのは, Julia側ではアルゴリズムの最適化をしっかりやっているからである.
#
# コンパイラによる最適化ではなく, 数学的なアルゴリズムの最適化なので, 原理的には同様の最適化は任意のプログラミング言語において可能である.
#
# だから, そういう意味では, 「Juliaが速い」のではなく, 「採用したアルゴリズムが速い」のである.

# %% [markdown]
# ### 最適化の中身の紹介
#
# しかし, アルゴリズムの最適化は「有理函数の分子と分母をHorner法でべた書きする」というステップを含んでいる.  最大で分子分母がそれぞれ31次と32次の多項式の場合に分子分母をHorner法でべた書きしたコードを書く必要がある. (Taylor展開で計算する部分については最大で37次の多項式をHorner法でべた書きする必要がある.)
#
# そのようなコードを人間が書くことは通常ないし, 人間に書かせることはバグの原因にもなるのでやってはいけないことである.
#
# MITでの宿題の答えではJuliaがLispのような完全なマクロ機能を持っていることを使って, そのようなべた書きのコードを自動的に生成している.

# %% [markdown]
# ### コードの自動生成機能の重要性
#
# 最適化の過程で必要になる試行錯誤の作業においてはそのようなコードの自動生成機能は必須である.
#
# なぜならば, Float64またはComplex{Float64}の場合に, 各複素数ごとに連分数展開をどれだけ長く取れば十分な精度が得られるかを確認し, 領域を適切に区切って, 領域ごとに使用する連分数展開の長さを決めるという作業が必要になるからだ.
#
# Horner法でべた書きしたコードと定義通りに連分数展開しているコードでは誤差の出方が異なる. だから, 十分な精度が得られる次数を決定するためには, Horner法でべた書きしたケースを直接調べる必要がある.
#
# そのような作業においてHorner法によるべた書きのコードを自動生成できないならば, 手間の分量は大変なことになってしまうだろう.

# %% [markdown]
# ### 教訓
#
# このように原理的に任意のプログラミング言語で可能な最適化であっても, 各プログラミング言語で使える環境の良し悪しによって, 「人間が最適化をやる気になれるかどうか」は大きく変わってしまうのだ.
#
# おそらくそのような理由でFortranで書かれたライブラリは十分に最適化されなかったのだろう.
#
# 以上のストーリーの教訓はこうである. 
#
# * コンパイラの最適化性能だけが, 現実に使用されているライブラリの実行速度を決めるのではない.
# * ライブラリを書いた人がシャープな最適化をやる気になれるという条件が非常に重要である.
# * そのためには, 数式処理系の支援が得られたり, 視覚化を駆使した試行錯誤が容易だったり, コードの自動生成が容易であったりすることが重要になる.
#
# Juliaはそれらの条件をすべて満たしていたおかげで最適化の作業が容易になり, 宿題の模範解答に収まる内容で, 既存のライブラリよりも数倍速いコードを実現することができたのである.

# %% [markdown]
# ### **「怠惰な人間」という巨大因子**を無視してはいけない
#
# 実践的な状況における実効的な速さについて考えるときには, コンパイラの最適化の性能が主な因子であると思っては**いけない**. 
#
# **「怠惰な人間」という巨大因子**を無視している「速さ」に関する議論はナンセンスである. 
#
# コンパイラの最適化の性能のみに特化したベンチマークテストを行っただけでは, 実践的な状況における実効的な速さについては大したことは分からない.

# %% [markdown]
# ### uncorrelated氏による連分数のベンチマークテスト
#
# 以上のような話を私はツイッターでよくしていた.
#
# 最近(2020年9月28日), 私がよく話題にしていた[MITの宿題の答えの一部分をコピーしたと思われるコードを含むベンチマークテスト](https://github.com/uncorrelated/ExpIntCF)を発表している人(uncorrelated氏)を発見した. その内容を見たところ, 以上で述べたような考え方を無視しているだけではなく, やっていることがあまりにも滅茶苦茶であった. そのような代物が**まともなもの扱い**されることは有害だと思ったので, このノートを書くことにした.

# %% [markdown]
# ### 補足と感謝
#
# このノートを書いた真の狙いは, MITの宿題の答え及びその帰結について詳しく説明することによって, 複雑なコードの自動生成可能性がどれだけ重要であるかをみんなに知ってもらうことである.  MITの宿題の模範解答のコードはuncorrelated氏のコードよりも100倍から1000倍も速い. 
#
# そのような桁違いに遅いコードでベンチマークテストを行ったuncorrelated氏のおかげで, それなりに有益なノートを書くことができたことには感謝の意を示したい. どうもありがとう!

# %% [markdown]
# ## MITでの講義の宿題の解答のコード
#
# 次のセルのコードは本質的に
#
# * https://nbviewer.jupyter.org/github/stevengj/18S096/blob/iap2017/pset3/pset3-solutions.ipynb
#
# からのコピーである.  ただし, Julia v1 以上に対応するために
#
# * [Polynomials.jl/src/polynomials/Poly.jl](https://github.com/JuliaMath/Polynomials.jl/blob/6def89946e386085349310b581f8286c2093855f/src/polynomials/Poly.jl)
# * https://github.com/JuliaLang/julia/blob/master/HISTORY.md#breaking-changes
#
# を参考にしながら修正してある:
#
# * `eulergamma` の代わりに `MathConstants.eulergamma` を使うようにした.
# * `Polynomials.Poly` の代わりに `Polynomials.PolyCompat.Poly` を使用するようにした.
# * `@macrocall` の文法の breaking change に対応した. (`Base.LineNumberNode(@__LINE__, Symbol(@__FILE__)),` を挿入した.)
# * 一ヶ所, `@E₁_cf` となっている部分を `@E₁_cf64` に修正した.
#
# 函数 $E_1(z)$ の実装結果の定義域は $\real(z) > 0$ である.
#
# 函数 $E_1(z)$ は $\real(z)^2 + 0.233\imag(z)^2 \ge 7.84 = 2.8^2$ ならば連分数展開で計算され, それ以外の場合にはTaylor展開で計算される.  連分数展開は「無限遠近傍」で効率的な数値計算法であり, 原点近傍では効率が悪いことに注意せよ.

# %%
# https://github.com/mitmath/18S096/blob/iap2017/pset3/pset3-solutions.ipynb

# SOLUTION code
# n coefficients of the Taylor series of E₁(z) + log(z), in type T:
function E₁_taylor_coefficients(::Type{T}, n::Integer) where T<:Number
    n < 0 && throw(ArgumentError("$n ≥ 0 is required"))
    n == 0 && return T[]
    n == 1 && return T[-MathConstants.eulergamma]
    # iteratively compute the terms in the series, starting with k=1
    term::T = 1
    terms = T[-MathConstants.eulergamma, term]
    for k=2:n
        term = -term * (k-1) / (k * k)
        push!(terms, term)
    end
    return terms
end

# SOLUTION code
macro E₁_taylor64(z, n::Integer)
    c = E₁_taylor_coefficients(Float64, n)
    taylor = Expr(:macrocall, Symbol("@evalpoly"),
        Base.LineNumberNode(@__LINE__, Symbol(@__FILE__)),
        :t, c...)
    quote
        let t = $(esc(z))
            $taylor - log(t)
        end
    end
end

# compute E₁ via n terms of the continued-fraction expansion, implemented
# in the simplest way:
function E₁_cf(z::Number, n::Integer)
    # starting with z seems to give many fewer terms for intermediate |z| ~ 3
    cf::typeof(inv(z)) = z
    for i = n:-1:1
        cf = z + (1+i)/cf
        cf = 1 + i/cf
    end
    return exp(-z) / (z + inv(cf))
end

# SOLUTION code
# for numeric-literal coefficients: simplify to a ratio of two polynomials:
import Polynomials
# return (p,q): the polynomials p(x) / q(x) corresponding to E₁_cf(x, a...),
# but without the exp(-x) term
function E₁_cfpoly(n::Integer, ::Type{T}=BigInt) where T<:Real
    q = Polynomials.PolyCompat.Poly(T[1])
    p = x = Polynomials.PolyCompat.Poly(T[0,1])
    for i = n:-1:1
        p, q = x*p+(1+i)*q, p # from cf = x + (1+i)/cf = x + (1+i)*q/p
        p, q = p + i*q, p     # from cf = 1 + i/cf = 1 + i*q/p
    end
    # do final 1/(x + inv(cf)) = 1/(x + q/p) = p/(x*p + q)
    return p, x*p + q
end
macro E₁_cf64(x, n::Integer)
    p,q = E₁_cfpoly(n, BigInt)
    evalpoly = Symbol("@evalpoly")
    num_expr = Expr(:macrocall, evalpoly, 
        Base.LineNumberNode(@__LINE__, Symbol(@__FILE__)),
        :t, Float64.(Polynomials.coeffs(p))...)
    den_expr = Expr(:macrocall, evalpoly,
        Base.LineNumberNode(@__LINE__, Symbol(@__FILE__)),
        :t, Float64.(Polynomials.coeffs(q))...)
    quote
        let t = $(esc(x))
            exp(-t) * $num_expr / $den_expr
        end
    end
end

# SOLUTION:
function E₁(z::Union{Float64,Complex{Float64}})
    x² = real(z)^2
    y² = imag(z)^2
    if x² + 0.233*y² ≥ 7.84 # use cf expansion, ≤ 30 terms
        if (x² ≥ 546121) & (real(z) > 0) # underflow
            return zero(z)
        elseif x² + 0.401*y² ≥ 58.0 # ≤ 15 terms
            if x² + 0.649*y² ≥ 540.0 # ≤ 8 terms
                x² + y² ≥ 4e4 && return @E₁_cf64 z 4
                return @E₁_cf64 z 8
            end
            return @E₁_cf64 z 15
        end
        return @E₁_cf64 z 30
    else # use Taylor expansion, ≤ 37 terms
        r² = x² + y²
        return r² ≤ 0.36 ? (r² ≤ 2.8e-3 ? (r² ≤ 2e-7 ? @E₁_taylor64(z,4) :
                                                       @E₁_taylor64(z,8)) :
                                         @E₁_taylor64(z,15)) :
                          @E₁_taylor64(z,37)
    end
end
E₁(z::Union{T,Complex{T},Rational{T},Complex{Rational{T}}}) where T<:Integer = E₁(float(z))

# %% [markdown]
# ### 連分数で計算する部分の解説
#
# 函数 $E_1(z)$ の定義域は $\real(z) > 0$ であり, 函数 $E_1(z)$ は $\real(z)^2 + 0.233\imag(z)^2 \ge 7.84 = 2.8^2$ ならば連分数展開で計算され, それ以外の場合にはTaylor展開で計算される. 
#
# 使用される連分数の定義は `E₁_cf(z)` を見ればわかる.  コードでは見難いので数式で表示させてみよう. 

# %%
@syms z
E₁_cf(z, 4)

# %% [markdown]
# これは $e^{-z}$ の因子を除けば有理函数である. その分子と分母は次のようになる.

# %%
E₁_cf(z, 4).simplify()

# %% [markdown]
# 連分数を定義通りにループを回して計算すると遅くなってしまう. 以下のように有理函数の分子分母をHorner法のコードをべた書きすると速くなる.  この例ではそれをJuliaのマクロ(コードの自動生成)を使って実現している.

# %%
@E₁_cf64(z, 4)

# %% [markdown]
# 以下の場合には分子が31次, 分母が32次の多項式になっている(連分数展開).

# %%
p, q = E₁_cfpoly(30, BigInt)
display(p)
display(q)

# %%
@vars z
print("p(z) = ")
show(@evalpoly z Float64.(Polynomials.coeffs(p))...)
print("\n\nq(z) = ")
show(@evalpoly z Float64.(Polynomials.coeffs(q))...)

# %% [markdown]
# 以下は37次の多項式である(Taylor展開).

# %%
r = Polynomials.PolyCompat.Poly(E₁_taylor_coefficients(Float64, 37))

# %%
@vars z
print("r(z) = ")
show(@evalpoly z Float64.(Polynomials.coeffs(r))...)

# %% [markdown]
# 計算したいのは, Float64の場合なので, 適当に分割された領域ごとにどれだけ長い連分数(もしくはTaylor展開)を計算すれば十分な精度が得られるかを前もって計算して確認することができる.  上に引用したMITでの宿題の模範解答のコード `E₁(z)` はそのようにして作られている.

# %% [markdown]
# ## uncorrelated氏のコード
#
# 以下に引用したuncorrelated氏のコードは
#
# * https://github.com/uncorrelated/ExpIntCF/blob/master/e1_cf.jl
#
# より.  ただし, 函数名に `_uc` を付け加えてある.

# %%
# https://github.com/uncorrelated/ExpIntCF/blob/master/e1_cf.jl

function E₁_cf_uc(z::Number, n::Integer)
    cf::typeof(z) = z
    for i = n:-1:1
        cf = z + (1+i)/cf
        cf = 1 + i/cf
    end
    return exp(-z) / (z + inv(cf))
end

function E₁_cf_uc(z::Number, reltol=1e-12)
    for n = 1:1000
        s = E₁_cf_uc(z, n)
        d = E₁_cf_uc(z, 2n)
        if abs(s - d) <= reltol*abs(d)
            return d
        end
    end
    error("iteration limit exceeded!")
end

# %% [markdown]
# 以下に再引用するMITの宿題の模範解答の一部分(in the simple wayで実装されているので, このコードは速くない)と上に引用したuncorrelated氏の `E₁_cf_uc(z::Number, n::Integer)` はほぼ同じであることに注意せよ.  これが単なる偶然なわけがない.

# %%
# 再掲
# https://nbviewer.jupyter.org/github/stevengj/18S096/blob/iap2017/pset3/pset3-solutions.ipynb
# より

# compute E₁ via n terms of the continued-fraction expansion, implemented
# in the simplest way:
function E₁_cf(z::Number, n::Integer)
    # starting with z seems to give many fewer terms for intermediate |z| ~ 3
    cf::typeof(inv(z)) = z
    for i = n:-1:1
        cf = z + (1+i)/cf
        cf = 1 + i/cf
    end
    return exp(-z) / (z + inv(cf))
end

# %% [markdown]
# ### 最初のバージョンでは reltol=1e-16 になっていた！(笑)
#
# なんと驚くべきことに, それの最初のバージョン
#
# * https://github.com/uncorrelated/ExpIntCF/blob/031ab9925ad08870737e771c062c12df90ab9750/e1_cf.jl
#
# では, 
#
# ```julia
# function E₁_cf(z::Number, reltol=1e-16)
#     略
# end
# ```
#
# となっていた. これだと reltol が machine epsilon (Float64のeps())より小さくなってしまう(笑).

# %%
eps()

# %%
1e-16 < eps()

# %% [markdown]
# ### 連分数展開が効率的でない領域で連分数のコードをテストしている！
#
# さらに, 現在のバージョン
#
# * https://github.com/uncorrelated/ExpIntCF/blob/90669379f0f1ca8478ddc3677f78bccfbcebcbdd/e1_cf.jl
#
# においても,
#
# ```julia
# function MakeMatrix(s::Number, e::Number, length::Integer)
# 	x = range(s, e, length=length)
# 	return [E₁_cf(x+y*im) for y in x, x in x]
# end
#
# m = @time MakeMatrix(0.1, 5, 100)
# ```
#
# となっている.  これは無限遠の近くでの展開である連分数展開を原点に近い $0.1\le \real(z), \imag(z)\le 5$ で使用していることになる. そのような計算は恐ろしく効率が悪い.

# %% [markdown]
# ### 実践的には決して使われることがない遅いアルゴリズムで比較している！
#
# さらに, テストに使用されたアルゴリズムは遅いことが前もって分かっているアルゴリズムである. 
#
# 次のセルの函数 `E₁_cf_uc(z::Number, n::Integer)` はMITでの宿題の解答の函数 `E₁_cf(z::Number, n::Integer)` のほぼコピーである. なぜか完全なコピーになっておらず, `cf::typeof(inv(z)) = z` が `cf::typeof(z) = z` になっている.
#
# その函数はそこでは最適化前のコード扱いされている. `E₁_cf(z::Number, n::Integer)` は有限連分数を定義通りに計算する函数であり, 本質的に有限連分数で定義された有理函数の値を計算するコードになっている.  有利函数は多項式分の多項式であり, 多項式の数値計算はHorner法のコードをべた書きすると速くなる.  MITの宿題の答えではそれをJuliaのマクロを使って実現している.
#
# 実践的な状況での実効的な速さを比較したければ, そのマクロを使って実現されているコードのアルゴリズムで比較するべきだろう. 
#
# しかし, 問題はJuliaのように完全なマクロを持たないプログラミング言語では, Horner法のべた書き部分を実現する手間が増える. だから, その増えた手間の分も評価しないと不公平な比較になってしまう.
#
# ところが, uncorrelated氏はそのような実践的な状況における実効的な速さのフェアな比較をしようとせずに, わざわざ最初から遅いと分かっているコードで比較している. (後で具体的に確認するが, 確認されたケースでは, MITでの宿題の模範解答と比較すると100倍から1000倍もuncorrelated氏が採用したコードは遅い.)

# %% [markdown]
# **補足:** `cf::typeof(inv(z)) = z` が `cf::typeof(z) = z` に変更してしまったせいで, uncorrelated氏の `E₁_cf_uc(z::Number, n::Integer)` では `z` が整数型の場合の値を計算できない.  (MITでの宿題の模範解答の `E₁_cf(z::Number, n::Integer)` では計算できる.)

# %%
E₁_cf_uc(3, 4) # by uncorrelated

# %%
E₁_cf(3, 4) # solution of MIT homework

# %% [markdown]
# ## uncorrelated氏の計算法がどれだけ遅くなっているか
#
# uncorrelated氏の計算法がどれだけ遅くなっているかを確認してみよう.

# %% [markdown]
# ### 素朴な比較

# %%
# 次の函数は
# https://github.com/uncorrelated/ExpIntCF/blob/master/e1_cf.jl
# より. 函数名に `_uc` を付け加えてある
function MakeMatrix_uc(s::Number, e::Number, length::Integer)
    x = range(s, e, length=length)
    return [E₁_cf_uc(x+y*im) for y in x, x in x]
end

M_uc = MakeMatrix_uc(0.1, 5, 100)
print("uncorrelated's code:")
@btime MakeMatrix_uc(0.1, 5, 100);

# %%
function MakeMatrix_mit(s::Number, e::Number, length::Integer)
    x = range(s, e, length=length)
    return [E₁(x+y*im) for y in x, x in x]
end

M_mit = MakeMatrix_mit(0.1, 5, 100)
print("MIT homework solution:")
@btime MakeMatrix_mit(0.1, 5, 100);

# %% [markdown]
# 数値が似ていること, 特にメモリ割当の数値が完全に同じであることによって, 後者に `m` の文字が付け加わっていることを見落とさないようにせよ！(笑) 
#
# uncorrelated氏のコードとMITの宿題の模範解答のコードでは計算速度の__桁が3桁も違う！__ __秒とミリ秒の違いがある！__

# %% [markdown]
# 以下を見ればわかるように, uncorrelated氏のコードとMITの宿題の答えの計算結果の差の相対誤差は1e-13未満になっている. uncorrelated氏はreltol=1e-12と設定していたので, これは妥当な結果だろう. 
#
# ヒートマップのカラーバーのメモリは「2つの値の比と1の差の絶対値の常用対数」である. 色が明るい場所ほど違いが大きい.

# %%
max_relerr = maximum(@. abs(M_uc/M_mit - 1))
@show max_relerr
x = y = range(0.1, 5, length=100)
heatmap(x, y, log10.(abs.(M_uc./M_mit .- 1)); size=(500, 450))

# %% [markdown]
# ### 連分数展開が有効な領域での比較
#
# 上の比較において, MITの宿題の解答の E₁_cf64(z) は原点の近傍ではTaylor展開を使って効率的に計算するようになっている.  一方, uncorrelated氏の計算法では原点の近傍でも効率の悪い連分数展開を使っている.  これだと余りにもuncorrelated氏にとって不利過ぎるので, 実部が $3$ 以上の確実に連分数展開が有効な領域で速度を比較してみよう.

# %%
N_uc = MakeMatrix_uc(3, 20, 100)
print("uncorrelated's code:")
@btime MakeMatrix_uc(3, 20, 100);

# %%
N_mit = MakeMatrix_mit(3, 20, 100)
print("MIT homework solution:")
@btime MakeMatrix_mit(3, 20, 100);

# %% [markdown]
# 連分数展開が確実に有効な領域においても, 計算速度が__2桁も違う！__

# %%
max_relerr = maximum(@. abs(N_uc/N_mit - 1))
@show max_relerr
x = y = range(0.1, 5, length=100)
heatmap(x, y, log10.(abs.(N_uc./N_mit .- 1)); size=(500, 450))

# %% [markdown]
# ヒートマップのカラーバーのメモリは「2つの値の比と1の差の絶対値の常用対数」である.　色が明るい場所ほど違いが大きい.

# %% [markdown]
# ### 結論: 論外！
#
# 以上のようにuncorrelated氏によるベンチマークテストは論外な内容であると言ってよい.

# %% [markdown]
# 実践的に使用可能なレベルのコードよりも百倍から千倍も遅いコードでベンチマークテストをやることにどういう意味があるのだろうか?
#
# 本当にやるべき作業はJuliaで実現したアルゴリズムと同じアルゴリズムをCまたはFortranで実装して, Juliaよりどれだけ速くなるかを確認することだと思われる. 私は実際に速くなると信じている. そのようなコードは貴重なライブラリとして広く配布する価値があるかもしれない.
#
# しかし, 適切なアルゴリズムの調整には試行錯誤が必須であり, そのような試行錯誤は現代的な「何でも揃っているプログラミング言語環境」でないと現実にはやる気になれないだろう. 
#
# さらに, このノートに引用したアルゴリズムを手間暇かけてCやFortranに移植しても, 「何でも揃っている側のプログラミング言語環境」におけるさらなるアルゴリズムの改良によって時代遅れになってしまうかもしれない. 
#
# CまたはFortranに移植する過程では, CまたはFortranで直接書くには面倒な部分のコードを「何でも揃っている側のプログラミング言語環境」で自動生成するようにした方がよいだろう. 
#
# そういうことができなければ, CまたはFortranでHorner法のべた書きが必要になる. 誰かがんばって試してくれ！(笑)

# %%
