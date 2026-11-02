# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
using Statistics
using Random

function verify_coin_paradox(n_coins = 1_000_000)
    frequencies = Float64[]
    
    for _ in 1:n_coins
        # 4回のコイントス (true = 表, false = 裏)
        flips = rand(Bool, 4)
        
        denom = 0 # 表(Heads)の直後にトスが行われた回数
        num = 0   # 表(Heads)の直後のトスが表だった回数
        
        # 4回目の直後にはトスがないため、1〜3回目をチェック
        for i in 1:3
            if flips[i]
                denom += 1
                if flips[i+1]
                    num += 1
                end
            end
        end
        
        # 少なくとも1回は「表の後のトス」が発生した場合のみ記録
        # (4回連続裏などの場合は相対頻度が計算できないため除外)
        if denom > 0
            push!(frequencies, num / denom)
        end
    end
    
    # 記録された相対頻度の平均を計算
    return mean(frequencies)
end

# 再現性のためのシード設定
Random.seed!(42)

# シミュレーション実行
result = verify_coin_paradox()
println("100万回のシミュレーションによる平均相対頻度: ", round(result, digits=5))

# %%
