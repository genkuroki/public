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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# # 5次式の2次式と3次式の積への因数分解の仕方
#
# * 黒木玄

# %%
using SymPy

# Override the Base.show definition of SymPy.jl:
# https://github.com/JuliaPy/SymPy.jl/blob/29c5bfd1d10ac53014fa7fef468bc8deccadc2fc/src/types.jl#L87-L105

@eval SymPy function Base.show(io::IO, ::MIME"text/latex", x::SymbolicObject)
    print(io, as_markdown("\\displaystyle " * sympy.latex(x, mode="plain", fold_short_frac=false)))
end
@eval SymPy function Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{Sym})
    function toeqnarray(x::Vector{Sym})
        a = join(["\\displaystyle " * sympy.latex(x[i]) for i in 1:length(x)], "\\\\")
        """\\left[ \\begin{array}{r}$a\\end{array} \\right]"""
    end
    function toeqnarray(x::AbstractArray{Sym,2})
        sz = size(x)
        a = join([join("\\displaystyle " .* map(sympy.latex, x[i,:]), "&") for i in 1:sz[1]], "\\\\")
        "\\left[ \\begin{array}{" * repeat("r",sz[2]) * "}" * a * "\\end{array}\\right]"
    end
    print(io, as_markdown(toeqnarray(x)))
end

# %% [markdown]
# ## 一般論

# %%
@vars a b c d e p q r s t x

# %% [markdown]
# $t\ne 0$ と仮定する.

# %%
f = x^5 + p*x^4 + q*x^3 + r*x^2 + s*x + t

# %%
g = x^2 + a*x + b

# %%
h = x^3 + c*x^2 + d*x + e

# %% [markdown]
# $f$ を $g$ で割ったときの商と余りを求める.

# %%
quotient, remainder = sympy.div(f, g, x)
sympy.poly.([quotient, remainder], x)

# %% [markdown]
# $h$ が商に一致することと $c,d,e$ が次のようになることは同値.

# %%
C, D, E = [quotient.coeff(x, k) for k in 2:-1:0]

# %% [markdown]
# $be = t\ne 0$ ならば $b\ne 0$ でなければいけない.
#
# そのとき, 余りが0になることと以下の2つが0になることは同値である.

# %%
eqrem0 = remainder.coeff(x, 0)/b |> expand

# %%
eqrem1 = remainder.coeff(x, 1)

# %% [markdown]
# 与えられた5次でモニックな整数係数多項式 $f$ の整数係数の2次式 $g$ と3次式 $h$ の積への分解を求めるためには, まず $be=t$ を満たす整数の組 $(b, e)$ の各々について, `eqrem0` を0にする整数 $a$ が存在するかどうかを確認し, 存在する場合に `eqrem1` が0になるかどうかも確認し, `eqrem0` と `eqrem1` の両方を $0$ にするような整数の組 $(a,b)$ を見付ければよい.  (見付からない場合には整数係数の2次式と整数係数の3次式の積に分解不可能なことが証明されたことになる.)

# %% [markdown]
# ## 例1
#
# $g = x^2 + 3x - 2$, $h = (x+1)^3-2 = x^3+3x^2+3x-1$ の場合にどうなるか.

# %%
GG = x^2 + 3x - 2

# %%
HH = (x+1)^3-2 |> expand

# %%
FF = GG*HH

# %%
FF = expand(FF)

# %% [markdown]
# 既存の函数を使えば因数分解は易しい.

# %%
factor(FF)

# %% [markdown]
# 係数を $P,Q,R,S,T$ に代入する.

# %%
_, P, Q, R, S, T = [FF.coeff(x, k) for k in 5:-1:0]

# %% [markdown]
# 余りの定数項が0になるという条件.

# %%
Eqrem0 = eqrem0(p=>P, q=>Q, r=>R, s=>S, t=>T)

# %% [markdown]
# 余りの1次の項が0になるという条件.

# %%
Eqrem1 = eqrem1(p=>P, q=>Q, r=>R, s=>S, t=>T)

# %% [markdown]
# $be = 2$ となる整数の組 $(b, e)$ 全体について, 最下段の方程式が整数解を持つ場合を探す.

# %%
[Eqrem0(b=>k) |> factor for k in (1, 2, -1, -2)]

# %% [markdown]
# 最後の $b=-2$ の場合にのみ整数解 $a=3$ が存在することがわかった.

# %%
BB, AA = -2, 3

# %% [markdown]
# この場合が実際に解になっていることの確認.

# %%
C, D, E

# %%
CC, DD, EE = C(a=>AA, b=>BB, p=>P), D(a=>AA, b=>BB, p=>P, q=>Q), E(a=>AA, b=>BB, p=>P, q=>Q, r=>R, t=>T)

# %%
Eqrem1(a=>AA, b=>BB)

# %%
sol = (x^2 + AA*x + BB)*(x^3 + CC*x^2 + DD*x + EE)

# %%
sol - FF |> simplify

# %% [markdown]
# ## 例2
#
# $f = x^5 - x^4 - 1$ の因数分解.

# %%
FF = x^5 - x^4 - 1

# %% [markdown]
# 既存の函数を使えば易しい.

# %%
factor(FF)

# %% [markdown]
# よくある解答例では天下り的に「$x^3$ を足して引く」という方法を使っている.
#
# $$
# \begin{aligned}
# x^5 - x^4 - 1
# &=
# x^5 - x^4 + x^3 - x^3 - 1
# \\ &=
# (x^5 - x^4 + x^3) - (x^3 + 1)
# \\ &=
# x^3(x^2-x+1) -(x+1)(x^2-x+1)
# \\ &=
# (x^2-x+1)(x^3-x-1).
# \end{aligned}
# $$
#
# 以下の方法であればそういう天下り的な操作が無用になる.

# %% [markdown]
# 係数を $P,Q,R,S,T$ に代入する.

# %%
_, P, Q, R, S, T = [FF.coeff(x, k) for k in 5:-1:0]

# %% [markdown]
# 余りの定数項が0になるという条件.

# %%
Eqrem0 = eqrem0(p=>P, q=>Q, r=>R, s=>S, t=>T)

# %% [markdown]
# 余りの1次の項が0になるという条件.

# %%
Eqrem1 = eqrem1(p=>P, q=>Q, r=>R, s=>S, t=>T)

# %% [markdown]
# $be=-1$ となる整数の組 $(b,d)$ の各々について, 最下段の方程式が整数解を持つかどうかを確認する.

# %%
[(B = k; E = T/k; Eqrem0(b=>B, e=>E) |> factor) for k in (1, -1)]

# %% [markdown]
# $b=\pm 1$ の両方の場合に整数解 $a=-1$ を持つことがわかった.
#
# これらが下から2段目の方程式も満たしているかを確認する.

# %%
BB, AA = 1, -1
Eqrem1(a=>AA, b=>BB)

# %%
BB, AA = -1, -1
Eqrem1(a=>AA, b=>BB)

# %% [markdown]
# $(a, b) = -1, 1$ の場合のみが解になっていることがわかった.

# %% [markdown]
# 以下は検算である.

# %%
AA, BB = -1, 1
CC, DD, EE = C(a=>AA, b=>BB, p=>P), D(a=>AA, b=>BB, p=>P, q=>Q), E(a=>AA, b=>BB, p=>P, q=>Q, r=>R, t=>T)

# %%
sol = (x^2 + AA*x + BB)*(x^3 + CC*x^2 + DD*x + EE)

# %%
sol - FF |> simplify

# %% [markdown]
# ## 例3
#
# $f = x^5 + 16x + 32$ の因数分解.

# %%
FF = x^5 + 16x + 32

# %% [markdown]
# 既存の函数を使えば易しい.

# %%
factor(FF)

# %% [markdown]
# 係数を $P,Q,R,S,T$ に代入する.

# %%
_, P, Q, R, S, T = [FF.coeff(x, k) for k in 5:-1:0]

# %% [markdown]
# 余りの定数項が0になるという条件.

# %%
Eqrem0 = eqrem0(p=>P, q=>Q, r=>R, s=>S, t=>T)

# %% [markdown]
# 余りの1次の項が0になるという条件.

# %%
Eqrem1 = eqrem1(p=>P, q=>Q, r=>R, s=>S, t=>T)

# %% [markdown]
# $b$ として $32$ の約数全体を動かして, `eqrem0` を0にする整数 $a$ が存在する場合を探す.

# %%
[[k, Eqrem0(b=>k) |> factor] for k in ((2 .^ (0:5))..., (-2 .^ (0:5))...)] |> (x -> stack(x; dims=1))

# %% [markdown]
# 方程式 `eqrem0 = 0` は $b=4,-2$ のときにのみ整数解 $a=2$ を持つことがわかった.

# %%
BB, AA = 4, 2
Eqrem1(a=>AA, b=>BB)

# %%
BB, AA = -2, 2
Eqrem1(a=>AA, b=>BB)

# %% [markdown]
# 方程式 `eqrem0 = eqrem1 = 0` の整数解 $a = 2$, $b=4$ が見付かった.

# %%
AA, BB = 2, 4
CC, DD, EE = C(a=>AA, b=>BB, p=>P), D(a=>AA, b=>BB, p=>P, q=>Q), E(a=>AA, b=>BB, p=>P, q=>Q, r=>R, t=>T)

# %%
sol = (x^2 + AA*x + BB)*(x^3 + CC*x^2 + DD*x + EE)

# %%
sol - FF |> simplify

# %% [markdown]
# ## 例4
#
# $f = x^5 + 3x^4 - 2x^3 - 2x^2 - 6x + 4$ の因数分解

# %%
FF = x^5 + 3x^4 - 2x^3 - 2x^2 - 6x + 4

# %%
factor(FF)

# %% [markdown]
# $f$ が $x^k + a$ 型の因子を持つと決め打ちできるなら, そのように $f$ を整理すれば容易に因数分解できる.
#
# $$
# \begin{aligned}
# x^5 + 3x^4 - 2x^3 - 2x^2 - 6x + 4
# &=
# (x^5 + 3x^4 - 2x^3) - (2x^2 + 6x - 4)
# \\ &=
# x^3(x^2 + 3x - 2) - 2(x^2 + 3x - 2)
# \\ &=
# (x^3-2)(x^2 + 3x - 2).
# \end{aligned}
# $$

# %%
factor(FF)

# %%
_, P, Q, R, S, T = [FF.coeff(x, k) for k in 5:-1:0]

# %%
eqrem0

# %%
Eqrem0 = eqrem0(p=>P, q=>Q, r=>R, s=>S, t=>T)

# %%
eqrem1

# %%
Eqrem1 = eqrem1(p=>P, q=>Q, r=>R, s=>S, t=>T)

# %%
[Eqrem0(b=>k, e=>T/k) |> factor for k in (1, 2, 4, -1, -2, -4)]

# %%
BB, AA = -2, 3
Eqrem1(a=>AA, b=>BB)

# %%
AA, BB = 3, -2
CC, DD, EE = C(a=>AA, b=>BB, p=>P), D(a=>AA, b=>BB, p=>P, q=>Q), E(a=>AA, b=>BB, p=>P, q=>Q, r=>R, t=>T)

# %%
sol = (x^2 + AA*x + BB)*(x^3 + CC*x^2 + DD*x + EE)

# %%
sol - FF |> simplify

# %%
