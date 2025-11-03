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
xstep = 0.0001
x = 0:xstep:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

# 図の作成
fig = Figure(size=(600, 400))

# サブプロット2: 複数の破線パターン
ax2 = Axis(fig[1,1], title="Multiple Dash Patterns CairoMakie x=$x")
lines!(ax2, x, y1, label="dash", linewidth=2, linestyle=:dash)
lines!(ax2, x, y2, label="dash", linewidth=2, linestyle=:dash)
lines!(ax2, x, y3, label="dash", linewidth=2, linestyle=:dash)
axislegend(ax2)
xlims!(ax2, 0, 10)

# PDF形式で保存（この際にApple製品での表示問題が発生する可能性あり）
save("dash_plot2_CairoMakie_$(xstep).pdf", fig)

# PNG形式でも保存（こちらは問題なく表示される）
save("dash_plot2_CairoMakie_$(xstep).png", fig)

fig

# %%
using Plots: Plots
Plots.gr(titlefontsize=12)

xstep = 0.0001
x = 0:xstep:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

P = Plots.plot(title="Multiple Dash Patterns by Plots gr() x=$x")
Plots.plot!(x, y1; ls=:dash, label="dash")
Plots.plot!(x, y2; ls=:dash, label="dash")
Plots.plot!(x, y3; ls=:dash, label="dash")
Plots.savefig("dash_plot2_Plots_gr_$xstep.png")
Plots.savefig("dash_plot2_Plots_gr_$xstep.pdf")
P

# %%
using Plots: Plots
Plots.pgfplotsx(titlefontsize=14)

# x = 0:0.0001:10 とすると以下の処理には非常に時間がかかってしまう。
# その意味でも x = 0:0.0001:10 とするのは止めた方が良い。
xstep = 0.01
x = 0:xstep:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

P = Plots.plot(title="Multiple Dash Patterns by Plots pgfplotsx() x=$x")
Plots.plot!(x, y1; ls=:dash, label="dash")
Plots.plot!(x, y2; ls=:dash, label="dash")
Plots.plot!(x, y3; ls=:dash, label="dash")
Plots.plot!(tex_output_standalone=true)
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep)_.tex")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).png")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).pdf")
P

# %%
using Plots: Plots
Plots.pgfplotsx(titlefontsize=14)

# x = 0:0.0001:10 とすると以下の処理には非常に時間がかかってしまう。
# その意味でも x = 0:0.0001:10 とするのは止めた方が良い。
xstep = 0.005
x = 0:xstep:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

P = Plots.plot(title="Multiple Dash Patterns by Plots pgfplotsx() x=$x")
Plots.plot!(x, y1; ls=:dash, label="dash")
Plots.plot!(x, y2; ls=:dash, label="dash")
Plots.plot!(x, y3; ls=:dash, label="dash")
Plots.plot!(tex_output_standalone=true)
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep)_.tex")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).png")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).pdf")
P

# %%
using Plots: Plots
Plots.pgfplotsx(titlefontsize=14)

# x = 0:0.0001:10 とすると以下の処理には非常に時間がかかってしまう。
# その意味でも x = 0:0.0001:10 とするのは止めた方が良い。
xstep = 0.001
x = 0:xstep:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

P = Plots.plot(title="Multiple Dash Patterns by Plots pgfplotsx() x=$x")
Plots.plot!(x, y1; ls=:dash, label="dash")
Plots.plot!(x, y2; ls=:dash, label="dash")
Plots.plot!(x, y3; ls=:dash, label="dash")
Plots.plot!(tex_output_standalone=true)
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep)_.tex")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).png")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).pdf")
P

# %%
using Plots: Plots
Plots.pgfplotsx(titlefontsize=14)

# x = 0:0.0001:10 とすると以下の処理には非常に時間がかかってしまう。
# その意味でも x = 0:0.0001:10 とするのは止めた方が良い。
xstep = 0.0005
x = 0:xstep:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

P = Plots.plot(title="Multiple Dash Patterns by Plots pgfplotsx() x=$x")
Plots.plot!(x, y1; ls=:dash, label="dash")
Plots.plot!(x, y2; ls=:dash, label="dash")
Plots.plot!(x, y3; ls=:dash, label="dash")
Plots.plot!(tex_output_standalone=true)
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep)_.tex")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).png")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).pdf")
P

# %%
using Plots: Plots
Plots.pgfplotsx(titlefontsize=14)

# x = 0:0.0001:10 とすると以下の処理には非常に時間がかかってしまう。
# その意味でも x = 0:0.0001:10 とするのは止めた方が良い。
xstep = 0.0001
x = 0:xstep:10
y1 = sin.(x)
y2 = sin.(x .+ 0.2)
y3 = sin.(x .+ 0.4)

P = Plots.plot(title="Multiple Dash Patterns by Plots pgfplotsx() x=$x")
Plots.plot!(x, y1; ls=:dash, label="dash")
Plots.plot!(x, y2; ls=:dash, label="dash")
Plots.plot!(x, y3; ls=:dash, label="dash")
Plots.plot!(tex_output_standalone=true)
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep)_.tex")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).png")
Plots.savefig("dash_plot2_Plots_pgfplotsx_$(xstep).pdf")
P

# %%
