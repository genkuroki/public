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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using BenchmarkTools

# %%
"""
    nextcombination!(n, t, c = typeof(t)[min(t-1, i) for i in 1:t])

`[1,2,…,n]` からの重複無しの `t` 個の組み合わせ `c` をすべて生成したい.

`nextcombination!(n, t, c)` は配列で表現された組み合わせ `c` をその次の組み合わせに書き換えて, `c` を返す.

初期条件を `c = typeof(t)[min(t-1, i) for i in 1:t]` にすると, `binomial(n, t)` 回の `nextcombination!(n, t, c)` ですべての組み合わせが生成される.
"""
function nextcombination!(n, t, c = typeof(t)[min(t-1, i) for i in 1:t])
    t == 0 && return c
    @inbounds for i in t:-1:1
        c[i] += 1
        c[i] > (n - (t - i)) && continue
        for j in i+1:t
            c[j] = c[j-1] + 1
        end
        break
    end
    c
end

"""
    mycombinations!(n::Integer, t, c)

事前に割り当てられた組み合わせを格納する配列 `c` を使って, `[1,2,…,n]` からの重複無しの `t` 個の組み合わせのすべてを生成する生成子を返す.
"""
function mycombinations!(n::Integer, t, c)
    for i in 1:t c[i] = min(t - 1, i) end
    (nextcombination!(n, t, c) for _ in 1:binomial(n, t))
end

"""
    mycombinations!(a, t, c)

事前に割り当てられた組み合わせを格納する配列 `c` を使って, 配列 `a` からのインデックスに重複がない `t` 個の組み合わせのすべてを生成する生成子を返す.
"""
function mycombinations!(a, t, c)
    t < 0 && (t = length(a) + 1)
    (view(a, indices) for indices in mycombinations!(length(a), t, c))
end

"""
    mycombinations(x, t)

`x` が整数ならば `[1,2,…,x]` からの, `x` が配列ならば `x` からのインデックスに重複がない `t` 個の組み合わせのすべてを生成する生成子を返す.
"""
mycombinations(x, t) = mycombinations!(x, t, Vector{typeof(t)}(undef, t))

# %%
@doc nextcombination!

# %%
@doc mycombinations!

# %%
@doc mycombinations

# %%
# nextcombination!(n, t, c) は組み合わせを表す配列 c を次の組み合わせに書き変える.
# nextcombination!(n, t, c) の結果は最初の binomial(n, t) 個だけが有効.
n, t = 5, 3
@eval @show binomial($n, $t)
@show c = typeof(t)[min(t-1, i) for i in 1:t]
for i in 1:binomial(n, t)+2
    @eval @show $i, nextcombination!($n, $t, c)
end

# %%
using StatsBase

"""
    mytiedrank(X,
        rank = similar(X, Float64),
        place = similar(X, Int))

配列 `X` 内の値の順位(下から何番目に大きいか)を計算する.

値に重複がある場合には同順位にするために順位の平均を順位として採用する.

`rank`, `place` は内部で使う配列で, 事前に割り当てておけば, この函数はメモリ割り当てゼロで実行される.  `rank` が返り値になる.
"""
function mytiedrank(X,
        rank = similar(X, Float64),
        place = similar(X, Int))
    # Xの要素が整数ならば sortperm(X) を使った方がずっと速い.
    # そうなる理由は sortperm は sort! と違って sortperm_int_range も使用するから.
    # その代わりに sortperm だとメモリ割り当てが発生する.
    # この辺は Base 自体に改善の余地が残っているように思われる.
    #place = sortperm(X)
    place .= axes(place, 1)
    sort!(place; by = i -> @inbounds(X[i]))
    StatsBase._tiedrank!(rank, X, place)
end

# %%
@doc mytiedrank

# %%
using StatsBase

X = Float64[8, 5, 1, 2, 3, 3, 5, 5]
rank = similar(X, Float64)
place = similar(X, Int)

@show X
a = @show mytiedrank(X, rank, place)
b = @show mytiedrank(X)
c = @show tiedrank(X)
a == b == c

# %%
X = rand(1.0:100, 100000)
rank = similar(X, Float64)
place = similar(X, Int)
a = @btime mytiedrank(X, rank, place)
b = @btime tiedrank(X)
a == b

# %%
"""
    exact_mann_whitney_u(x, y,
        xy = Vector{promote_type(eltype(x), eltype(y))}(undef, length(x)+length(y)),
        rankxy = similar(xy, Float64),
        place = similar(xy, Int),
        c = similar(x, Int, min(length(x), length(y)))
    )

データ `x, y` に関するMann-WhitneyのU検定のexact版を素朴なアルゴリズム(遅い)で実行し, `(u, u_less, u_greater, p_twosided, p_less, p_greater)` を返す.  ここで,

* `u` は `x`, `y` の `length` が小さい方のU統計量達.
* `u_less`, `u_greater` はそれぞれ `u`, `length(x)*length(y) - u` の小さい方と大きい方.
* `p_twosided` は両側検定のP値.
* `p_less` は `u_less` 以下になる確率.
* `p_greater` は `u_greater` 以上になる確率.

`p_twosided` は `p_less + p_greater` と1の小さい方になる.

`xy`, `rankxy`, `place`, `comb` は内部で使用する配列である. それらを事前に割り当てておけば, この函数はメモリ割り当てゼロで実行される.

この函数の実装では素朴にすべての組み合わせに関する和を計算する方法を採用しているので, `x`, `y` のサイズが大きくなると極めて遅くなる.  より洗練されたアルゴリズムについては以下の文献を参照せよ.

* https://pure.tue.nl/ws/files/1330412/200012964.pdf
* https://www.sciencedirect.com/science/article/pii/S1672022916000358

注意: この函数が返す両側検定のP値はR言語の [exactRankTests::wilcox.exact](https://rdrr.io/cran/exactRankTests/src/R/wilcox.exact.R) や [coin::wilcox_test](https://rdrr.io/cran/coin/src/R/LocationTests.R) と同じ値になる.  それらは素朴な実装を採用しているこの函数よりも計算が速い.

注意: この函数が返す両側検定のP値は [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) の `pvalue(ExactMannWhitneyUTest(x, y))` とは一般に異なる.  そうなる理由は HypothesisTests.jl が `[x; y]` 内の値に重複(tie)がある場合に分布の左右対称性が崩れることに配慮していないからである.  さらに現時点の HypothesisTests.jl (v0.10.10) の `pvalue(ExactMannWhitneyUTest(x, y))` ではメモリ割り当てがかなり発生する.

注意: この函数が返す両側検定のP値は [scipy.stats.mannwhitneyu](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html) とは一般に異なる.  その理由は仕様的に "No correction is made for ties." となっているからである.
"""
function exact_mann_whitney_u(x, y,
        xy = Vector{promote_type(eltype(x), eltype(y))}(undef, length(x)+length(y)),
        rankxy = similar(xy, Float64),
        place = similar(xy, Int),
        comb = similar(x, Int, min(length(x), length(y)))
    )
    # Assume length(x) ≤ length(y)
    if length(x) > length(y)
        return exact_mann_whitney_u(y, x, xy, rankxy, place, comb)
    end
    
    # Initialization
    m, n = length(x), length(y)
    xy[1:m] .= x
    xy[m+1:m+n] .= y

    # Calculation of ranks
    rankxy = mytiedrank(xy, rankxy, place)
    
    # Calculation of P-values
    r = sum(@view rankxy[1:m]) # rank sum of the data x,y
    r_less, r_greater = minmax(r, m*(m+n+1) - r)
    num_less, num_greater = 0, 0
    for c in mycombinations!(rankxy, m, comb)
        R = sum(c) # rank sum of the combination c
        num_less    += R ≤ r_less
        num_greater += R ≥ r_greater
    end
    p_less, p_greater = (num_less, num_greater) ./ binomial(m+n, m)
    p_twosided = min(1, p_less + p_greater)
    
    # Calculation of U-statistics
    u, u_less, u_greater = (r, r_less, r_greater) .- m*(m+1)÷2
    
    # Returns
    u, u_less, u_greater, p_twosided, p_less, p_greater
end

# %%
@doc exact_mann_whitney_u

# %%
using BenchmarkTools
using HypothesisTests

# https://www.jstor.org/stable/2533173
x = Float64[1,2,1,1,1,1,1,1,1,1,2,4,1,1]
y = Float64[3,3,4,3,1,2,3,1,1,5,4]

xy = Vector{promote_type(eltype(x), eltype(y))}(undef, length(x)+length(y))
rankxy = similar(xy, Float64)
place = similar(xy, Int)
comb = similar(x, Int, min(length(x), length(y)))

a = @btime exact_mann_whitney_u($x, $y)
b = @btime exact_mann_whitney_u($x, $y, $xy, $rankxy, $place, $comb)
c = @btime pvalue(ExactMannWhitneyUTest($x, $y))
@show a b c;

# %%
using RCall
@rimport exactRankTests as exactRankTests
using SciPy

# %%
x, y = Float64[1, 2], Float64[0, 5, 3, 4]
@show x y
@show m, n = length(x), length(y)
@show exact_mann_whitney_u(x, y)
@show exact_mann_whitney_u(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y))
@show pvalue(ExactMannWhitneyUTest(y, x))
@show rcopy(exactRankTests.wilcox_exact(x, y))[:p_value]
@show rcopy(exactRankTests.wilcox_exact(y, x))[:p_value]
@show SciPy.stats.mannwhitneyu(x, y; method="exact")
@show SciPy.stats.mannwhitneyu(y, x; method="exact");

# %%
x, y = Float64[1, 2], Float64[2, 2, 3, 4]
@show x y
@show m, n = length(x), length(y)
@show exact_mann_whitney_u(x, y)
@show exact_mann_whitney_u(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y))
@show pvalue(ExactMannWhitneyUTest(y, x))
@show rcopy(exactRankTests.wilcox_exact(x, y))[:p_value]
@show rcopy(exactRankTests.wilcox_exact(y, x))[:p_value]
@show SciPy.stats.mannwhitneyu(x, y; method="exact")
@show SciPy.stats.mannwhitneyu(y, x; method="exact");

# %%
x, y = Float64[1, 2], Float64[2, 2, 2, 4, 5]
@show x y
@show m, n = length(x), length(y)
@show exact_mann_whitney_u(x, y)
@show exact_mann_whitney_u(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y))
@show pvalue(ExactMannWhitneyUTest(y, x))
@show rcopy(exactRankTests.wilcox_exact(x, y))[:p_value]
@show rcopy(exactRankTests.wilcox_exact(y, x))[:p_value]
@show SciPy.stats.mannwhitneyu(x, y; method="exact")
@show SciPy.stats.mannwhitneyu(y, x; method="exact");

# %% tags=[]
x, y = Float64[1, 2, 2], Float64[2, 2, 2, 4, 5, 6, 7]
@show x y
@show m, n = length(x), length(y)
@show exact_mann_whitney_u(x, y)
@show exact_mann_whitney_u(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y))
@show pvalue(ExactMannWhitneyUTest(y, x))
@show rcopy(exactRankTests.wilcox_exact(x, y))[:p_value]
@show rcopy(exactRankTests.wilcox_exact(y, x))[:p_value]
@show SciPy.stats.mannwhitneyu(x, y; method="exact")
@show SciPy.stats.mannwhitneyu(y, x; method="exact");

# %%
x = Float64[11, 15,  9,  4, 34, 17, 18, 14, 12, 13, 26, 31]
y = Float64[34, 31, 35, 29, 28, 12, 18, 30, 14, 22, 10]
@show x y
@show m, n = length(x), length(y)
@show exact_mann_whitney_u(x, y)
@show exact_mann_whitney_u(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y))
@show pvalue(ExactMannWhitneyUTest(y, x))
@show rcopy(exactRankTests.wilcox_exact(x, y))[:p_value]
@show rcopy(exactRankTests.wilcox_exact(y, x))[:p_value]
@show SciPy.stats.mannwhitneyu(x, y; method="exact")
@show SciPy.stats.mannwhitneyu(y, x; method="exact");

# %%
x = Float64[ 4,  7,  8,  9, 13, 13, 17, 11] |> sort
y = Float64[23,  6,  3, 24, 17, 14, 24, 29, 13, 33] |> sort
@show x y
@show m, n = length(x), length(y)
@show exact_mann_whitney_u(x, y)
@show exact_mann_whitney_u(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y))
@show pvalue(ExactMannWhitneyUTest(y, x))
@show rcopy(exactRankTests.wilcox_exact(x, y))[:p_value]
@show rcopy(exactRankTests.wilcox_exact(y, x))[:p_value]
@show SciPy.stats.mannwhitneyu(x, y; method="exact")
@show SciPy.stats.mannwhitneyu(y, x; method="exact");

# %%
x = Float64[ 4,  7,  8,  9, 13, 13, 17, 11] |> sort
y = Float64[23,  6,  3, 24, 17, 14, 24, 29, 13, 33] |> sort
@show x y
@show m, n = length(x), length(y)
println()

xy = Vector{promote_type(eltype(x), eltype(y))}(undef, length(x)+length(y))
rankxy = similar(xy, Float64)
place = similar(xy, Int)
comb = similar(x, Int, min(length(x), length(y)))

a = @btime exact_mann_whitney_u($x, $y, $xy, $rankxy, $place, $comb)
b = @btime exact_mann_whitney_u($x, $y)
c = @btime exact_mann_whitney_u($y, $x, $xy, $rankxy, $place, $comb)
d = @btime exact_mann_whitney_u($y, $x)
e = @btime pvalue(ExactMannWhitneyUTest($x, $y))
f = @btime rcopy(exactRankTests.wilcox_exact($x, $y))[:p_value]
g = @btime SciPy.stats.mannwhitneyu($x, $y; method="exact")
@show a b c d e f g;

# %%
x = Float64[4, 7, 8, 9, 11, 12, 13, 17]
y = Float64[3, 6, 10, 14, 16, 23, 24, 25, 29, 33]
@show x y
@show m, n = length(x), length(y)
println()

xy = Vector{promote_type(eltype(x), eltype(y))}(undef, length(x)+length(y))
rankxy = similar(xy, Float64)
place = similar(xy, Int)
comb = similar(x, Int, min(length(x), length(y)))

a = @btime exact_mann_whitney_u($x, $y, $xy, $rankxy, $place, $comb)
b = @btime exact_mann_whitney_u($x, $y)
c = @btime exact_mann_whitney_u($y, $x, $xy, $rankxy, $place, $comb)
d = @btime exact_mann_whitney_u($y, $x)
e = @btime pvalue(ExactMannWhitneyUTest($x, $y))
f = @btime rcopy(exactRankTests.wilcox_exact($x, $y))[:p_value]
g = @btime SciPy.stats.mannwhitneyu($x, $y; method="exact")
@show a b c d e f g;

# %%
