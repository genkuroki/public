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
function pvalue_exact_mann_whitney_u_test(x, y,
        xy = Vector{promote_type(eltype(x), eltype(y))}(undef, length(x)+length(y)),
        rankxy = similar(xy, Float64),
        place = similar(xy, Int),
        c = similar(x, Int, min(length(x), length(y)))
    )
    # Initialization
    m, n = length(x), length(y)
    xy[1:m] .= x
    xy[m+1:m+n] .= y
    N = m + n
    place .= 1:N

    # Calculation of ranks
    sort!(place; by = i->xy[i])
    i = 1
    @inbounds while i ≤ N
        j = i
        vi = xy[place[i]]
        while (j + 1 ≤ N) && (vi == xy[place[j + 1]])
            j += 1
        end
        if j > i
            t = j - i + 1
            rk = sum(i:j) / t
            for k in i:j
                rankxy[place[k]] = rk
            end
        else
            rankxy[place[i]] = i
        end
        i = j + 1
    end
    
    # Calculation of the two-sided exact P-value
    l = min(m, n)
    r = l == m ? sum(@view rankxy[1:m]) : sum(@view rankxy[m+1:m+n])
    r_le, r_ge = minmax(r, m*n + l*(l+1) - r)
    le, ge = 0, 0
    for comb in mycombinations!(rankxy, l, c)
        R = sum(comb)
        le += R ≤ r_le
        ge += R ≥ r_ge
    end
    min(1, (le + ge)/binomial(N, l))
end

# %%
using BenchmarkTools
using HypothesisTests
using RCall
R"library(coin)"
R"library(exactRankTests)"
using DataFrames
using SciPy

# %%
x, y = [1, 2], [0, 5, 3, 4]
m, n = length(x), length(y)
@show pvalue_exact_mann_whitney_u_test(x, y)
@show pvalue_exact_mann_whitney_u_test(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y))

# %%
@rput x y
R"""
wilcox.exact(x, y)
"""

# %%
score = [x; y]
gender = [fill(1, length(x)); fill(0, length(y))]
df = DataFrame(score=score, gender=gender)
@show df
@rput df
R"""
df$gender <- as.factor(df$gender)
wilcox_test(score ~ gender, distribution = "exact", data=df)
"""

# %%
SciPy.stats.mannwhitneyu(x, y; method="exact")

# %%
x, y = [1, 2], [2, 2, 3, 4]
m, n = length(x), length(y)
@show pvalue_exact_mann_whitney_u_test(x, y)
@show pvalue_exact_mann_whitney_u_test(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y));

# %%
@rput x y
R"""
wilcox.exact(x, y)
"""

# %%
score = [x; y]
gender = [fill(1, length(x)); fill(0, length(y))]
df = DataFrame(score=score, gender=gender)
@show df
@rput df
R"""
df$gender <- as.factor(df$gender)
wilcox_test(score ~ gender, distribution = "exact", data=df)
"""

# %%
SciPy.stats.mannwhitneyu(x, y; method="exact")

# %%
x, y = [1, 2], [2, 2, 2, 4, 5]
m, n = length(x), length(y)
@show pvalue_exact_mann_whitney_u_test(x, y)
@show pvalue_exact_mann_whitney_u_test(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y));

# %%
@rput x y
R"""
wilcox.exact(x, y)
"""

# %%
SciPy.stats.mannwhitneyu(x, y; method="exact")

# %% tags=[]
x, y = [1, 2, 2], [2, 2, 2, 4, 5, 6, 7]
m, n = length(x), length(y)
@show pvalue_exact_mann_whitney_u_test(x, y)
@show pvalue_exact_mann_whitney_u_test(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y));

# %%
@rput x y
R"""
wilcox.exact(x, y)
"""

# %%
score = [x; y]
gender = [fill(1, length(x)); fill(0, length(y))]
df = DataFrame(score=score, gender=gender)
@show df
@rput df
R"""
df$gender <- as.factor(df$gender)
wilcox_test(score ~ gender, distribution = "exact", data=df)
"""

# %%
SciPy.stats.mannwhitneyu(x, y; method="exact")

# %%
x = [11, 15,  9,  4, 34, 17, 18, 14, 12, 13, 26, 31]
y = [34, 31, 35, 29, 28, 12, 18, 30, 14, 22, 10]
m, n = length(x), length(y)
@show pvalue_exact_mann_whitney_u_test(x, y)
@show pvalue_exact_mann_whitney_u_test(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y));

# %%
@rput x y
R"""
wilcox.exact(x, y)
"""

# %%
score = [x; y]
gender = [fill(1, length(x)); fill(0, length(y))]
df = DataFrame(score=score, gender=gender)
@show df
@rput df
R"""
df$gender <- as.factor(df$gender)
wilcox_test(score ~ gender, distribution = "exact", data=df)
"""

# %%
SciPy.stats.mannwhitneyu(x, y; method="exact")

# %%
@show x = [ 4,  7,  8,  9, 13, 13, 17, 11] |> sort
@show y = [23,  6,  3, 24, 17, 14, 24, 29, 13, 33] |> sort
m, n = length(x), length(y)
println()
@show pvalue_exact_mann_whitney_u_test(x, y)
@show pvalue_exact_mann_whitney_u_test(y, x)
@show pvalue(ExactMannWhitneyUTest(x, y));

# %%
@rput x y
R"""
wilcox.exact(x, y)
"""

# %%
score = [x; y]
gender = [fill(1, length(x)); fill(0, length(y))]
df = DataFrame(score=score, gender=gender)
@show df
@rput df
R"""
df$gender <- as.factor(df$gender)
wilcox_test(score ~ gender, distribution = "exact", data=df)
"""

# %%
SciPy.stats.mannwhitneyu(x, y; method="exact")

# %%
@show x = [ 4,  7,  8,  9, 13, 13, 17, 11] |> sort
@show y = [23,  6,  3, 24, 17, 14, 24, 29, 13, 33] |> sort
@show m, n = length(x), length(y)
println()

xy = Vector{promote_type(eltype(x), eltype(y))}(undef, length(x)+length(y))
rankxy = similar(xy, Float64)
place = similar(xy, Int)
c = similar(x, Int, min(length(x), length(y)))

a = @btime pvalue_exact_mann_whitney_u_test($x, $y, $xy, $rankxy, $place, $c)
b = @btime pvalue_exact_mann_whitney_u_test($x, $y)
c = @btime pvalue_exact_mann_whitney_u_test($y, $x, $xy, $rankxy, $place, $c)
d = @btime pvalue_exact_mann_whitney_u_test($y, $x)
e = @btime pvalue(ExactMannWhitneyUTest($x, $y))
@show a b c d e;

# %%
