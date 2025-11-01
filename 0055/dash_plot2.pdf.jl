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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
using Makie, CairoMakie

# データの準備
x = 0:0.1:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

# 図の作成
fig = Figure(size=(600, 400))

# サブプロット2: 複数の破線パターン
ax2 = Axis(fig[1,1], title="Multiple Dash Patterns")
lines!(ax2, x, y1, label="dash", linewidth=2, linestyle=:dash)
lines!(ax2, x, y2, label="dash", linewidth=2, linestyle=:dash)
lines!(ax2, x, y3, label="dash", linewidth=2, linestyle=:dash)
axislegend(ax2)
xlims!(ax2, 0, 10)

# PDF形式で保存（この際にApple製品での表示問題が発生する可能性あり）
save("dash_plot.pdf", fig)

# PNG形式でも保存（こちらは問題なく表示される）
save("dash_plot.png", fig)

# %%
using Makie, CairoMakie

# データの準備
x = 0:0.0001:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

# 図の作成
fig = Figure(size=(600, 400))

# サブプロット2: 複数の破線パターン
ax2 = Axis(fig[1,1], title="Multiple Dash Patterns")
lines!(ax2, x, y1, label="dash", linewidth=2, linestyle=:dash)
lines!(ax2, x, y2, label="dash", linewidth=2, linestyle=:dash)
lines!(ax2, x, y3, label="dash", linewidth=2, linestyle=:dash)
axislegend(ax2)
xlims!(ax2, 0, 10)

# PDF形式で保存（この際にApple製品での表示問題が発生する可能性あり）
save("dash_plot2.pdf", fig)

# PNG形式でも保存（こちらは問題なく表示される）
save("dash_plot2.png", fig)

# %%
