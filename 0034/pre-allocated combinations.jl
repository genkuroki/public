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
using Combinatorics

function sumcombs(a, t)
    s = 0
    for c in combinations(a, t) # by Combinatorics.jl
        s += sum(c)
    end
    s
end

function mysumcombs(a, t, c = Vector{typeof(t)}(undef, t))
    s = 0
    for v in mycombinations!(a, t, c)
        s += sum(v)
    end
    s
end

a = rand(1:10, 20)
t = 10
c = Vector{typeof(t)}(undef, t)

@show sumcombs(a, t)
@show mysumcombs(a, t)
@show mysumcombs(a, t, c);

# %%
using BenchmarkTools

# %%
# Combinatorics.jl の combinations を使うとメモリ割当てが発生する.

@time sumcombs(a, t)
@time sumcombs(a, t)
@time sumcombs(a, t)
@btime sumcombs($a, $t)

# %%
# mycombinations! を使うとメモリ割当てを大幅に抑制できる.

@time mysumcombs(a, t)
@time mysumcombs(a, t)
@time mysumcombs(a, t)
@btime mysumcombs($a, $t)

# %%
# 事前割当てされた c と mycombinations! を使うとメモリ割当てをゼロにできる.

@time mysumcombs(a, t, c)
@time mysumcombs(a, t, c)
@time mysumcombs(a, t, c)
@btime mysumcombs($a, $t, $c)

# %%
@code_warntype mysumcombs(a, t, c)

# %%
