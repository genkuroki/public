# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.7.0-DEV
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# # 正弦函数の13次多項式による近似
#
# 区間 $[-\pi/4, \pi/4]$ で $\sin x$ を効率よく近似する13次多項式 $F_{13}(x)$ を計算する. 以下で $F_{10}(x)$ は10次多項式.
#
# * `Remez_abs`: $\sin(x) - F_{13}(x)$ のsupノルムを最小化
# * `Remez_rel`: $\dfrac{\sin(x) - F_{13}(x)}{\sin(x)}$ のsupノルムを最小化)
# * `Remez_rel1`: $\dfrac{\sin(x) - x - x^3 F_{10}(x)}{\sin(x)}$ のsupノルムを最小化)
# * `Remez_rel2`: $\dfrac{\sin(x) - x - x^3 F_{10}(x)}{x}$ のsupノルムを最小化)
#
# 最後の `Remez_rel2` の場合の $x + x^3 F_{10}(x)$ が [Julia の trig.jl](https://github.com/JuliaLang/julia/blob/master/base/special/trig.jl#L54) で使われている13次多項式に非常に近い.
#
# さらに, 低次の係数からFloat64での値に順次固定して求めて行くと, さらに近くなる.

# %% [markdown] toc=true
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Remez.jlの使用例" data-toc-modified-id="Remez.jlの使用例-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Remez.jlの使用例</a></span></li><li><span><a href="#低次の係数からFloat64での値に順次固定して求めてみる" data-toc-modified-id="低次の係数からFloat64での値に順次固定して求めてみる-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>低次の係数からFloat64での値に順次固定して求めてみる</a></span><ul class="toc-item"><li><span><a href="#weight-1/x-(minimaximize-the-pseudo-relative-error)" data-toc-modified-id="weight-1/x-(minimaximize-the-pseudo-relative-error)-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>weight 1/x (minimaximize the pseudo relative error)</a></span></li><li><span><a href="#weight-1/sin(x)-(minimaximize-the-exact-relative-error)" data-toc-modified-id="weight-1/sin(x)-(minimaximize-the-exact-relative-error)-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>weight 1/sin(x) (minimaximize the exact relative error)</a></span></li></ul></li><li><span><a href="#以下は汚い試行錯誤の残骸" data-toc-modified-id="以下は汚い試行錯誤の残骸-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>以下は汚い試行錯誤の残骸</a></span></li></ul></div>

# %% [markdown]
# https://github.com/simonbyrne/Remez.jl

# %%
using Remez
using Quadmath
using Polynomials
using Plots

relerr(x, y) = iszero(x - y) ? x - y : (x - y)/y 

# %%
?ratfn_minimax

# %% [markdown]
# ## Remez.jlの使用例

# %%
@time N0, D0, E0, X0 = ratfn_minimax(sin, (-big(π)/4, big(π)/4), 13, 0)
sleep(0.1)
c0 = @. Float128(round(N0; digits=50))
c0[begin+1:2:end] |> display
sin_Remez_abs = Polynomial(c0)

# %%
function chebyshev_approx(d, F::Polynomial{CT}) where CT
    P = F
    for k in degree(F):-1:d+1
        degree(P) < k && continue
        a = coeffs(P)[end]
        T = ChebyshevT([zeros(CT, k); a/CT(2)^(k-1)])
        P = P - T
    end
    P
end

function maclaurin_sinpio4(d, CT=BigFloat)
    s = [isodd(k) ? (-1)^((k-1)÷2)*(big(π)/4)^k/factorial(CT(k)) : zero(CT) for k in 0:d]
    Polynomial(s)
end

function maclaurin_sin(d, CT=BigFloat)
    s = [isodd(k) ? (-1)^((k-1)÷2)/factorial(CT(k)) : zero(CT) for k in 0:d]
    Polynomial(s)
end

function ordinaryscale(F::Polynomial{CT}) where CT
    c = coeffs(F) ./ [(CT(π)/4)^k for k in 0:degree(F)]
    Polynomial(c)
end

S1 = chebyshev_approx(13, maclaurin_sinpio4(61)) |> ordinaryscale
S2 = chebyshev_approx(13, maclaurin_sinpio4(63)) |> ordinaryscale
S3 = chebyshev_approx(13, maclaurin_sinpio4(65)) |> ordinaryscale
@show [
    S1 - S2
    S2 - S3
]
sleep(0.1)

s = Float128.(coeffs(S2))
s[begin+1:2:end] |> display
sin_genkuroki_abs = Polynomial(s)

# %%
coeffs(sin_genkuroki_abs - sin_Remez_abs)[begin+1:2:end]

# %%
coeffs(sin_Remez_abs)[begin+1:2:end] .|> Float64

# %%
@time N1, D1, E1, X1 = ratfn_minimax(sin, (-big(π)/4, big(π)/4), 13, 0, (x,y)->1/(eps(BigFloat)+abs(y)))
sleep(0.1)
c1 = @. Float128(round(N1; digits=50))
c1[begin+1:2:end] |> display
sin_Remez_rel = Polynomial(c1)

# %%
coeffs(sin_Remez_rel)[begin+1:2:end] .|> Float64

# %%
sinmxox3(x) = iszero(x) ? -1/big(6) : (sin(x) - x)/x^3 # sin(x) minus X over x^3
w_sinmxox3(x, y) = iszero(x) ? x : abs(x^3/sin(x))

@time N2, D2, E2, X2 = ratfn_minimax(sinmxox3, (-big(π)/4, big(π)/4), 10, 0, w_sinmxox3)
sleep(0.1)
c2 = @. Float128(round(N2; digits=50))
c2 = Float128[0; 1; 0; c2]
c2[begin+1:2:end] |> display
sin_Remez_rel1 = Polynomial(c2)

# %%
coeffs(sin_Remez_rel1)[begin+1:2:end] .|> Float64

# %% [markdown]
# * https://github.com/JuliaLang/julia/blob/master/base/special/trig.jl#L54
# * https://github.com/JuliaMath/openlibm/blob/master/src/k_sin.c#L51
# * https://github.com/freebsd/freebsd-src/blob/main/lib/msun/src/k_sin.c#L50

# %%
DS0 =  1.0
DS1 = -1.66666666666666324348e-01
DS2 =  8.33333333332248946124e-03
DS3 = -1.98412698298579493134e-04
DS4 =  2.75573137070700676789e-06
DS5 = -2.50507602534068634195e-08
DS6 =  1.58969099521155010221e-10

sin_julia = Polynomial(Float128[0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6])

# %%
coeffs(sin_Remez_rel1 - sin_julia)[begin+1:2:end]

# %%
[0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6][begin+1:2:end]

# %%
Float64.(coeffs(sin_Remez_rel1)[begin+1:2:end])

# %%
(Float64[0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6] - 
    Float64.(coeffs(sin_Remez_rel1)))[begin+1:2:end]

# %%
sinmxox3(x) = iszero(x) ? -1/big(6) : (sin(x) - x)/x^3 # sin(x) minus X over x^3

@time N3, D3, E3, X3 = ratfn_minimax(sinmxox3, (-big(π)/4, big(π)/4), 10, 0, (x,y)->x^2)
sleep(0.1)
c3 = @. Float128(round(N3; digits=50))
c3 = Float128[0; 1; 0; c3]
c3[begin+1:2:end] |> display
sin_Remez_rel2 = Polynomial(c3)

# %%
coeffs(sin_Remez_rel2 - sin_julia)[begin+1:2:end]

# %%
Float64.(coeffs(sin_Remez_rel1)[begin+1:2:end])

# %%
[0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6][begin+1:2:end]

# %%
([0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6] - 
    Float64.(coeffs(sin_Remez_rel2)))[begin+1:2:end]

# %%
x = range(eps(Float128), π/4; length=1000)
y_Remez_abs  = sin_Remez_abs.(x)
y_Remez_rel  = sin_Remez_rel.(x)
y_Remez_rel1 = sin_Remez_rel1.(x)
y_Remez_rel2 = sin_Remez_rel2.(x)
y_julia      = sin_julia.(x)
y            = sin.(x)

error_Remez_abs  = y_Remez_abs  - y .|> Float64
error_Remez_rel  = y_Remez_rel  - y .|> Float64
error_Remez_rel1 = y_Remez_rel1 - y .|> Float64
error_Remez_rel2 = y_Remez_rel2 - y .|> Float64
error_julia      = y_julia      - y .|> Float64

relerror_Remez_abs  = relerr.(y_Remez_abs , y) .|> Float64
relerror_Remez_rel  = relerr.(y_Remez_rel , y) .|> Float64
relerror_Remez_rel1 = relerr.(y_Remez_rel1, y) .|> Float64
relerror_Remez_rel2 = relerr.(y_Remez_rel2, y) .|> Float64
relerror_julia      = relerr.(y_julia     , y) .|> Float64

@show extrema(error_Remez_abs)
@show extrema(error_Remez_rel)
@show extrema(error_Remez_rel1)
@show extrema(error_Remez_rel2)
@show extrema(error_julia)
println()
@show extrema(relerror_Remez_abs)
@show extrema(relerror_Remez_rel)
@show extrema(relerror_Remez_rel1)
@show extrema(relerror_Remez_rel2)
@show extrema(relerror_julia)
;

# %%
plot(; legend=:bottom)
plot!(; title="10¹⁶ × absolute errors", titlefontsize=11)
plot!(x, 10^16*error_julia; label="Julia", ls=:solid)
plot!(x, 10^16*error_Remez_rel;  label="Remez_rel", ls=:dash)
plot!(x, 10^16*error_Remez_rel1; label="Remez_rel1", ls=:dot, lw=1.8)
plot!(x, 10^16*error_Remez_rel1; label="Remez_rel2", ls=:dash)
plot!(x, 10^16*error_Remez_abs;  label="Remez_abs", ls=:dashdot)

# %%
plot(; legend=:bottom)
plot!(; title="10¹⁶ × relative errors", titlefontsize=11)
plot!(x, 10^16*relerror_julia; label="Julia", ls=:solid)
plot!(x, 10^16*relerror_Remez_rel;  label="Remez_rel", ls=:dash)
plot!(x, 10^16*relerror_Remez_rel1; label="Remez_rel1", ls=:dot, lw=1.8)
plot!(x, 10^16*relerror_Remez_rel1; label="Remez_rel2", ls=:dash)
plot!(x, 10^16*relerror_Remez_abs;  label="Remez_abs", ls=:dashdot)

# %%
plot(; legend=:bottom)
plot!(; title="10¹⁶ × relative errors", titlefontsize=11)
plot!(x, 10^16*relerror_julia; label="Julia", ls=:solid)
plot!(x, 10^16*relerror_Remez_rel;  label="Remez_rel",  ls=:dash)
plot!(x, 10^16*relerror_Remez_rel1; label="Remez_rel1", ls=:dot, lw=1.8)
plot!(x, 10^16*relerror_Remez_rel1; label="Remez_rel2", ls=:dash)
plot!(x, 10^16*relerror_Remez_abs;  label="Remez_abs",  ls=:dashdot)
plot!(; ylim=(-0.08, 0.06))

# %%
plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors", titlefontsize=11)
plot!(x, 10^16*relerror_julia; label="Julia", ls=:solid)
plot!(x, 10^16*relerror_Remez_rel2; label="Remez_rel2", ls=:dash)

# %%
sin_Remez_rel2_F64 = convert(Polynomial{Float64}, sin_Remez_rel2)
sin_julia_F64 = Polynomial([0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6])

x = range(eps(Float128), π/4; length=1000)
y_Remez_rel2_F64 = sin_Remez_rel2_F64.(x)
y_julia_F64 = sin_julia_F64.(x)

relerror_Remez_rel2_F64 = relerr.(y_Remez_rel2_F64, y) .|> Float64
relerror_julia_F64 = relerr.(y_julia_F64, y) .|> Float64

plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors", titlefontsize=11)
plot!(x, 10^16*relerror_julia_F64; label="Julia", ls=:solid)
plot!(x, 10^16*relerror_Remez_rel2_F64; label="Remez_rel2_F64", ls=:dash)

# %% [markdown]
# ## 低次の係数からFloat64での値に順次固定して求めてみる

# %% [markdown]
# ### weight 1/x (minimaximize the pseudo relative error)

# %%
a1 = big"1.0"

# %%
r(x) = √abs(x)
F(x) = iszero(x) ? -1/typeof(x)(6) : (sin(x) - x*a1)/x^3
f(x) = F(r(x))

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 5, 0, (x,y)->abs(x))
sleep(0.1)

c3 = Float64.(N)
@show a3 = c3[1]
c3

# %%
F(x) = iszero(x) ? 1/big(120) : (sin(x) - x*(a1 + x^2*a3))/x^5
f(x) = F(r(x))

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 4, 0, (x,y)->abs(x)^2)
sleep(0.1)

c5 = Float64.(N)
@show a5 = c5[1]
c5

# %%
F(x) = iszero(x) ? -1/factorial(big(7)) : (sin(x) - x*(a1 + x^2*(a3 + x^2*a5)))/x^7
f(x) = F(r(x))

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 3, 0, (x,y)->abs(x)^3)
sleep(0.1)

c7 = Float64.(N)
@show a7 = c7[1]
c7

# %%
F(x) = iszero(x) ? 1/factorial(big(9)) : (sin(x) - x*(a1 + x^2*(a3 + x^2*(a5 + x^2*a7))))/x^9
f(x) = F(r(x))

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 2, 0, (x,y)->abs(x)^4)
sleep(0.1)

c9 = Float64.(N)
@show a9 = c9[1]
c9

# %%
F(x) = iszero(x) ? 1/factorial(big(11)) : (sin(x) - x*(a1 + x^2*(a3 + x^2*(a5 + x^2*(a7 + x^2*a9)))))/x^11
f(x) = F(r(x))

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 1, 0, (x,y)->abs(x)^5)
sleep(0.1)

c11 = Float64.(N)
@show a11 = c11[1]
c11

# %%
F(x) = iszero(x) ? 1/factorial(big(13)) : 
    (sin(x) - x*(a1 + x^2*(a3 + x^2*(a5 + x^2*(a7 + x^2*(a9 + x^2*a11))))))/x^13
f(x) = F(r(x))

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 0, 0, (x,y)->abs(x)^6)
sleep(0.1)

c13 = Float64.(N)
@show a13 = c13[1]
c13

# %%
a = Float64[0.0, a1, 0.0, a3, 0.0, a5, 0.0, a7, 0.0, a9, 0.0, a11, 0.0, a13]
a[begin+1:2:end]

# %%
a_julia = [0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6]
a_julia[begin+1:2:end]

# %%
(a - a_julia)[begin+1:2:end]

# %%
sin_a = Polynomial(a)
sin_j = Polynomial(a_julia)

x = range(eps(Float128), Float128(π)/4; length=1000)
y_a  = sin_a.(x)
y_j  = sin_j.(x)
y = sin.(x)

relerror_a = relerr.(y_a, y) .|> Float64
relerror_j = relerr.(y_j, y) .|> Float64
relerror_aj = relerr.(y_a, y_j) .|> Float64

@show extrema(relerror_a)
@show extrema(relerror_j)
@show extrema(relerror_aj)

plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors", titlefontsize=11)
plot!(x, 10^16*relerror_j; label="Julia", ls=:solid, lw=1.2)
plot!(x, 10^16*relerror_a; label="Mine",  ls=:dash, lw=1.2)

# %%
sin_a = Polynomial(a)
sin_j = Polynomial(a_julia)

x = range(eps(), π/4; length=1000)
y_a  = sin_a.(x)
y_j  = sin_j.(x)
y = sin.(Float128.(x))

relerror_a = relerr.(y_a, y) .|> Float64
relerror_j = relerr.(y_j, y) .|> Float64
relerror_aj = relerr.(y_a, y_j) .|> Float64

@show extrema(relerror_a)
@show extrema(relerror_j)
@show extrema(relerror_aj)

P = plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors (Julia)", titlefontsize=11)
plot!(x, 10^16*relerror_j; label="", lw=0.5)

Q = plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors (Mine)", titlefontsize=11)
plot!(x, 10^16*relerror_a; label="", lw=0.5)

plot(P, Q; size=(1000, 300))

# %% [markdown]
# ### weight 1/sin(x) (minimaximize the exact relative error)

# %%
a1 = big"1.0"

# %%
r(x) = √abs(x)
F(x) = iszero(x) ? -1/typeof(x)(6) : (sin(x) - x*a1)/x^3
f(x) = F(r(x))
W(x, y) = iszero(x) ? x : x^3/sin(x)
w(x, y) = W(r(x), y)

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 5, 0, w)
sleep(0.1)

@show N
c3 = Float64.(N)
@show a3 = c3[1]
c3

# %%
F(x) = iszero(x) ? -1/typeof(x)(6) : (sin(x) - x*(a1 + x^2*a3))/x^5
f(x) = F(r(x))
W(x, y) = iszero(x) ? x : x^5/sin(x)
w(x, y) = W(r(x), y)

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 4, 0, w)
sleep(0.1)

c5 = Float64.(N)
@show a5 = c5[1]
c5

# %%
F(x) = iszero(x) ? -1/factorial(big(7)) : (sin(x) - x*(a1 + x^2*(a3 + x^2*a5)))/x^7
f(x) = F(r(x))
W(x, y) = iszero(x) ? x : x^7/sin(x)
w(x, y) = W(r(x), y)

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 3, 0, w)
sleep(0.1)

c7 = Float64.(N)
@show a7 = c7[1]
c7

# %%
F(x) = iszero(x) ? 1/factorial(big(9)) : (sin(x) - x*(a1 + x^2*(a3 + x^2*(a5 + x^2*a7))))/x^9
f(x) = F(r(x))
W(x, y) = iszero(x) ? x : x^9/sin(x)
w(x, y) = W(r(x), y)

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 2, 0, w)
sleep(0.1)

c9 = Float64.(N)
@show a9 = c9[1]
c9

# %%
F(x) = iszero(x) ? 1/factorial(big(9)) : (sin(x) - x*(a1 + x^2*(a3 + x^2*(a5 + x^2*(a7 + x^2*a9)))))/x^11
f(x) = F(r(x))
W(x, y) = iszero(x) ? x : x^11/sin(x)
w(x, y) = W(r(x), y)

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 1, 0, w)
sleep(0.1)

c11 = Float64.(N)
@show a11 = c11[1]
c11

# %%
F(x) = iszero(x) ? 1/factorial(big(9)) : 
    (sin(x) - x*(a1 + x^2*(a3 + x^2*(a5 + x^2*(a7 + x^2*(a9 + x^2*a11))))))/x^13
f(x) = F(r(x))
W(x, y) = iszero(x) ? x : x^13/sin(x)
w(x, y) = W(r(x), y)

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 0, 0, w)
sleep(0.1)

c13 = Float64.(N)
@show a13 = c13[1]
c13

# %%
b = Float64[0.0, a1, 0.0, a3, 0.0, a5, 0.0, a7, 0.0, a9, 0.0, a11, 0.0, a13]
b[begin+1:2:end]

# %%
sin_b = Polynomial(b)

x = range(eps(Float128), Float128(π)/4; length=1000)
y_a  = sin_a.(x)
y_b  = sin_b.(x)
y_j  = sin_j.(x)
y = sin.(x)

relerror_a = relerr.(y_a, y) .|> Float64
relerror_b = relerr.(y_b, y) .|> Float64
relerror_j = relerr.(y_j, y) .|> Float64
relerror_aj = relerr.(y_a, y_j) .|> Float64
relerror_bj = relerr.(y_b, y_j) .|> Float64

@show extrema(relerror_a)
@show extrema(relerror_b)
@show extrema(relerror_j)
@show extrema(relerror_aj)
@show extrema(relerror_bj)

plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors", titlefontsize=11)
plot!(x, 10^16*relerror_j; label="Julia", ls=:solid, lw=1.2)
plot!(x, 10^16*relerror_a; label="My a",  ls=:dash, lw=1.2)
plot!(x, 10^16*relerror_b; label="My b",  ls=:dashdot, lw=1.2)

# %%
sin_a = Polynomial(a)
sin_b = Polynomial(b)
sin_j = Polynomial(a_julia)

x = range(eps(), π/4; length=1000)
y_a  = sin_a.(x)
y_b  = sin_b.(x)
y_j  = sin_j.(x)
y = sin.(Float128.(x))

relerror_a = relerr.(y_a, y) .|> Float64
relerror_b = relerr.(y_b, y) .|> Float64
relerror_j = relerr.(y_j, y) .|> Float64
relerror_aj = relerr.(y_a, y_j) .|> Float64
relerror_bj = relerr.(y_b, y_j) .|> Float64

@show extrema(relerror_a)
@show extrema(relerror_b)
@show extrema(relerror_j)
@show extrema(relerror_aj)
@show extrema(relerror_bj)

P = plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors (Julia)", titlefontsize=11)
plot!(x, 10^16*relerror_j; label="", lw=0.5)

Q = plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors (My b)", titlefontsize=11)
plot!(x, 10^16*relerror_b; label="", lw=0.5)

plot(P, Q; size=(1000, 300))

# %% [markdown]
# ## 以下は汚い試行錯誤の残骸

# %%
r(x) = √abs(x)
f(x) = iszero(x) ? -1/big(6) : (sin(r(x)) - r(x))/r(x)^3
w_f(x, y) = iszero(x) ? x : abs(r(x)^3/sin(r(x)))

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 5, 0, w_f)
sleep(0.1)

c = eltype(N3)[0; 1; 0; [isodd(k) ? N[(k+1)÷2] : 0 for k in 1:11]]
c = Float128.(c)
c[begin+1:2:end] |> display
sin_Remez_rel1′ = Polynomial(c)

# %%
sin_Remez_rel1′ - sin_Remez_rel1

# %%
r(x) = √abs(x)
f(x) = iszero(x) ? -1/big(6) : (sin(r(x)) - r(x))/r(x)^3
w_f(x, y) = iszero(x) ? x : abs(r(x)^3/sin(r(x)))

@time N, D, E, X = ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 5, 0, (x,y)->abs(x))
sleep(0.1)

c = eltype(N)[0; 1; 0; [isodd(k) ? N[(k+1)÷2] : 0 for k in 1:11]]
c = Float128.(c)
c[begin+1:2:end] |> display
sin_Remez_rel2′ = Polynomial(c)

# %%
sin_Remez_rel2′ - sin_Remez_rel2

# %%
using Remez: ratfn_leastsquares, ratfn_eval, find_extrema, winnow_extrema, ratfn_equal_deviation

function my_ratfn_minimax(f, interval, n, d,
                       w = (x,y)->BigFloat(1);
        epsbits = precision(BigFloat),
        threshold = 2^(-epsbits/3)
    )
    # We start off by finding a least-squares approximation. This
    # doesn't need to be perfect, but if we can get it reasonably good
    # then it'll save iterations in the refining stage.
    #
    # Least-squares approximations tend to look nicer in a minimax
    # sense if you evaluate the function at a big pile of Chebyshev
    # nodes rather than uniformly spaced points. These values will
    # also make a good grid to use for the initial search for error
    # extrema, so we'll keep them around for that reason too.

    # Construct the grid.
    lo = BigFloat(minimum(interval))
    hi = BigFloat(maximum(interval))

    local grid
    let
        mid = (hi+lo)/2
        halfwid = (hi-lo)/2
        nnodes = 16 * (n+d+1)
        grid = [ mid - halfwid * cospi(big(i)/big(nnodes)) for i=0:nnodes ]
    end

    # Find the initial least-squares approximation.
    (nc, dc) = ratfn_leastsquares(f, grid, n, d, w)

    # Threshold of convergence. We stop when the relative difference
    # between the min and max (winnowed) error extrema is less than
    # this.
    #
    # This is set to the cube root of machine epsilon on a more or
    # less empirical basis, because the rational-function case will
    # not converge reliably if you set it to only the square root.
    # (Repeatable by using the --test mode.) On the assumption that
    # input and output error in each iteration can be expected to be
    # related by a simple power law (because it'll just be down to how
    # many leading terms of a Taylor series are zero), the cube root
    # was the next thing to try.
    #epsbits = precision(BigFloat)
    #threshold = 2^(-epsbits/3)

    # Main loop.
    maxiters = 10^3
    iter = 0
    max_err = zero(BigFloat)
    extrema = Vector{Tuple{BigFloat, BigFloat}}()
    while iter < maxiters
        iter += 1
        iszero(mod(iter, 100)) && print(iter, " ")
        # Find all the error extrema we can.
        function compute_error(x)
            real_y = f(x)
            approx_y = ratfn_eval(nc, dc, x)
            return (approx_y - real_y) * w(x, real_y)
        end
        extrema = find_extrema(compute_error, grid)

        # Winnow the extrema down to the right number, and ensure they
        # have alternating sign.
        extrema = winnow_extrema(extrema, n+d+2)

        # See if we've finished.
        min_err = minimum([abs(y) for (x,y) = extrema])
        max_err = maximum([abs(y) for (x,y) = extrema])
        variation = (max_err - min_err) / max_err
        if variation < threshold
            return nc, dc, max_err, extrema
        end

        # If not, refine our function by equalising the error at the
        # extrema points, and go round again.
        (nc, dc) = ratfn_equal_deviation(f, map(x->x[1], extrema),
                                         n, d, max_err, w)
    end
    
    println("maxiters = $maxiters")
    return nc, dc, max_err, extrema
end

# %%
r(x) = √abs(x)
f(x) = iszero(x) ? -1/big(6) : (sin(r(x)) - r(x))/r(x)^3
w_f(x, y) = iszero(x) ? x : abs(r(x)^3/sin(r(x)))

setprecision(80)

@time N, D, E, X = my_ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 5, 0, w_f; 
    threshold = 2^(-precision(BigFloat)/7))
sleep(0.1)

c = eltype(N)[0; 1; 0; [isodd(k) ? N[(k+1)÷2] : 0 for k in 1:11]]
c[begin+1:2:end] |> display

sin_Remez = Polynomial(c)
sin_Remez |> display

setprecision(256);

# %%
r(x) = √abs(x)
f(x) = iszero(x) ? -1/big(6) : (sin(r(x)) - r(x))/r(x)^3
w_f(x, y) = iszero(x) ? x : abs(r(x)^3/sin(r(x)))

setprecision(64)

@time N, D, E, X = my_ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 5, 0, w_f; 
    threshold = 2^(-precision(BigFloat)/7))
sleep(0.1)

c = eltype(N)[0; 1; 0; [isodd(k) ? N[(k+1)÷2] : 0 for k in 1:11]]
c[begin+1:2:end] |> display

sin_Remez = Polynomial(c)
sin_Remez |> display

setprecision(256);

# %%
r(x) = √abs(x)
f(x) = iszero(x) ? -1/big(6) : (sin(r(x)) - r(x))/r(x)^3
w_f(x, y) = iszero(x) ? x : abs(r(x)^3/sin(r(x)))

setprecision(60)

@time N, D, E, X = my_ratfn_minimax(f, (big"0.0", (big(π)/4)^2), 5, 0, w_f; 
    threshold = 2^(-precision(BigFloat)/12))
sleep(0.1)

c = eltype(N)[0; 1; 0; [isodd(k) ? N[(k+1)÷2] : 0 for k in 1:11]]
c[begin+1:2:end] |> display
c[begin+1:2:end] .|> Float64 |> display

sin_Remez = Polynomial(c)
sin_Remez |> display

setprecision(256);

# %%
[0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6][begin+1:2:end] |> display

# %%
g1(x) = iszero(x) ? x : (sin(x) - x)/x
@time N, D, E, X = ratfn_minimax(g1, (-big(π)/4, big(π)/4), 12, 0)
sleep(0.1)
c = @. Float128(round(N; digits=50))
c = Float128[0; 1+c[1]; c[2:end]]
c[begin+1:2:end] |> display
sin_Remez_g1 = Polynomial(c)
sin_Remez_g1 |> display
y_g1 = sin_Remez_g1.(x)
relerror_Remez_g1 = relerr.(y_g1, y) .|> Float64

g2(x) = iszero(x) ? x : (sin(x) - x)/x^2
@time N, D, E, X = ratfn_minimax(g2, (-big(π)/4, big(π)/4), 11, 0, (x,y)->abs(x))
sleep(0.1)
c = @. Float128(round(N; digits=50))
c = Float128[0; 1; c]
c[begin+1:2:end] |> display
sin_Remez_g2 = Polynomial(c)
sin_Remez_g2 |> display
y_g2 = sin_Remez_g2.(x)
relerror_Remez_g2 = relerr.(y_g2, y) .|> Float64

g3(x) = iszero(x) ? -1/typeof(x)(6) : (sin(x) - x)/x^3
@time N, D, E, X = ratfn_minimax(g3, (-big(π)/4, big(π)/4), 10, 0, (x,y)->x^2)
sleep(0.1)
c = @. Float128(round(N; digits=50))
c = Float128[0; 1; 0; c]
c[begin+1:2:end] |> display
sin_Remez_g3 = Polynomial(c)
sin_Remez_g3 |> display
y_g3 = sin_Remez_g3.(x)
relerror_Remez_g3 = relerr.(y_g3, y) .|> Float64

plot(; legend=:bottomleft)
plot!(; title="10¹⁶ × relative errors (gₖ(x)=(sin(x)-x)/xᵏ)", titlefontsize=11)
plot!(x, 10^16*relerror_julia; label="Julia", ls=:solid)
plot!(x, 10^16*relerror_Remez_g1;  label="Remez_g1",  ls=:dash)
plot!(x, 10^16*relerror_Remez_g2;  label="Remez_g2",  ls=:dot, lw=1.8)
plot!(x, 10^16*relerror_Remez_g3;  label="Remez_g3",  ls=:dashdot)

# %%
cg2 = Float64.(coeffs(sin_Remez_g2))[begin+1:2:end]

# %%
cg3 = Float64.(coeffs(sin_Remez_g3))[begin+1:2:end]

# %%
cg2 == cg3

# %%
cjulia = [0.0, DS0, 0.0, DS1, 0.0, DS2, 0.0, DS3, 0.0, DS4, 0.0, DS5, 0.0, DS6][begin+1:2:end]

# %%
cg3 - cjulia

# %%
