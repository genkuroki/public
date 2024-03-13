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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # KP方程式のN-soliton解
#
# * 黒木玄
# * 2024-03-13
#
# See https://arxiv.org/abs/1208.2904

# %% [markdown]
# ## KP方程式の広田型N-soliton解の公式
#
# 与えるデータ: $u_i$, $v_i$, $c_i$ ($i=1,\ldots, N$)
#
# 独立変数: $t = (t_1, t_2, t_3, \ldots)$
#
# $\tau$ 函数:
#
# $$
# \tau(t) = \sum_{n=0}^N \sum_{1\le i_1<i_2<\cdots<i_n\le N}
# \prod_{1\le\mu<\nu\le n}\frac{(u_{i_\mu} - u_{i_\nu})(v_{i_\mu} - v_{i_\nu})}{(u_{i_\mu} - v_{i_\nu})(v_{i_\mu} - u_{i_\nu})}
# \prod_{\nu=1}^n c_{i_\nu} \exp\left(\sum_{m=1}^\infty (v_{i_\nu}^m - u_{i_\nu}^m) t_m\right).
# $$
#
# $n=0$ の項は $1$ であることに注意せよ.
#
# 簡単のため $t_1=x$, $t_2=y$, $t_3 = t$, $t_m=0$ ($m\ge 4$) とおく. そのとき,
#
# $$
# (D_x^4 + 3D_y^4 -4D_x D_t)\tau\bullet\tau = 0
# $$
#
# が成立している. ここで, $D$ 達は次のように定義される広田の$D$-operatorである:
#
# $$
# D_x^k f(x) \bullet g(x) =
# \left.\left(\frac{\partial}{\partial\varepsilon}\right)^k\,\right|_{\varepsilon=0} f(x+\varepsilon)g(x-\varepsilon) =
# \sum_{i=0}^k (-1)^{k-i}\binom{k}{i} f^{(i)}(x) g^{(k-i)}(x).
# $$
#
# このことから,
#
# $$
# u = 2\left(\frac{\partial}{\partial x}\right)^2 \log\tau
# $$
#
# とおくと, $u$ は次のKadomtsev-Petviashvili方程式を満たしていることが導かれる:
#
# $$
# \frac{\partial}{\partial x}\left(u_t - \frac{3}{2}uu_x - \frac{1}{4}u_{xxx}\right) - \frac{3}{4}u_{yy} = 0
# $$

# %%
using Combinatorics
using ForwardDiff
using Plots
default(fmt=:png, colorbar=false, tickfontsize=6)
using SymPy

function rescale_cgrad(x; f=x->x^2)
    cg = cgrad(x)
    cgrad(cg.colors, f.(cg.values))
end

nan2zero(x) = isnan(x) ? zero(x) : x

_eta(u, v, x, y, t) = (v - u) * x + (v^2 - u^2) * y + (v^3 - u^3) * t

function _soliton_exp_factor(u, v, c, i, x, y, t)
    n = length(i)
    prod_of_cs = prod(c[i[ν]] for ν in 1:n; init=1.0)
    sum_of_etas = sum(_eta(u[i[ν]], v[i[ν]], x, y, t) for ν in 1:n; init=0.0)
    prod_of_cs * exp(sum_of_etas)
end

function _soliton_tau_term(u, v, c, i, x, y, t)
    n = length(i)
    A = prod((u[i[μ]] - u[i[ν]]) * (v[i[μ]] - v[i[ν]]) / ((u[i[μ]] - v[i[ν]]) * (v[i[μ]] - u[i[ν]]))
            for μ in 1:n for ν in μ+1:n; init=1.0)
    A * _soliton_exp_factor(u, v, c, i, x, y, t)
end

function _soliton_tau(u, v, c, x, y, t)
    N = length(c)
    sum(_soliton_tau_term(u, v, c, i, x, y, t) for n in 0:N for i in combinations(1:N, n))
end

function make_soliton_tau(u, v, c)
    tau(x, y, t) = _soliton_tau(u, v, c, x, y, t)
    tau
end

make_soliton_tau(u, v) = make_soliton_tau(u, v, ones(length(u)))
make_soliton_tau_AB(k, P, c) = make_soliton_tau((P-k)/2, (P+k)/2, c)
make_soliton_tau_AB(k, P) = make_soliton_tau_AB(k, P, ones(length(k)))

first_derivative(f, x) = ForwardDiff.derivative(f, x)
second_derivative(f, x) = ForwardDiff.derivative(x -> ForwardDiff.derivative(f, x), x)

function tau2sol(tau)
    tau_x(x, y, t) = first_derivative(x -> tau(x, y, t), x)
    tau_xx(x, y, t) = second_derivative(x -> tau(x, y, t), x)
    sol(x, y, t) = 2(tau_xx(x, y, t) / tau(x, y, t) - (tau_x(x, y, t) / tau(x, y, t))^2)
    sol
end

function _make_soliton_solution(u, v, c)
    @syms x y t
    τ = _soliton_tau(u, v, c, x, y, t)
    τ_x = sympy.diff(τ, x)
    τ_xx = sympy.diff(τ_x, x)
    tau = eval(Meta.parse("(x, y, t) -> $τ"))
    tau_x = eval(Meta.parse("(x, y, t) -> $τ_x"))
    tau_xx = eval(Meta.parse("(x, y, t) -> $τ_xx"))
    sol(x, y, t) = 2(tau_xx(x, y, t) / tau(x, y, t) - (tau_x(x, y, t) / tau(x, y, t))^2)
    (; sol, tau, tau_x, tau_xx)
end

make_soliton_solution(u, v, c) = _make_soliton_solution(u, v, c).sol
make_soliton_solution(u, v) = make_soliton_solution(u, v, ones(length(u)))
make_soliton_solution_AB(k, P, c) = make_soliton_solution((P-k)/2, (P+k)/2, c)
make_soliton_solution_AB(k, P) = make_soliton_solution_AB(k, P, ones(length(k)))

function plot_sol_heatmap(sol, xs, ys, t; color=rescale_cgrad(:gist_earth))
    solmax = maximum(nan2zero(sol(x, y, t)) for x in xs, y in ys)
    heatmap(xs, ys, (x, y) -> sol(x, y, t); colorbar=false, color)
    plot!(; xguide="x", yguide="y")
    plot!(; clim=(-0.01solmax, solmax))
    plot!(; size=(300, 300))
end
    
function gif_sol_heatmap(sol, xs, ys, ts;
        fn="tmp.gif", fps=20, color=rescale_cgrad(:gist_earth)
    )
    solmax = maximum(nan2zero(sol(x, y, t)) for x in xs, y in ys, t in ts)
    anim = @animate for t in ts
        heatmap(xs, ys, (x, y) -> sol(x, y, t); colorbar=false, color)
        plot!(; xguide="x", yguide="y")
        plot!(; clim=(-0.01solmax, solmax))
        plot!(; size=(300, 300))
    end
    gif(anim, fn; fps)
end

function plot_sol_3d(sol, xs, ys, t; color=rescale_cgrad(:gist_earth))
    solmax = maximum(nan2zero(sol(x, y, t)) for x in xs, y in ys)
    surface(xs, ys, (x, y) -> sol(x, y, t); camera=(15, 80), colorbar=false, color)
    plot!(; xguide="x", yguide="y", zguide="u")
    plot!(; zlim=(-0.01solmax, solmax))
    plot!(; size=(400, 400), margin=-20Plots.mm)
end

function gif_sol_3d(sol, xs, ys, ts;
        fn="tmp.gif", fps=20, color=rescale_cgrad(:gist_earth)
    )
    solmax = maximum(nan2zero(sol(x, y, t)) for x in xs, y in ys, t in ts)
    anim = @animate for t in ts
        surface(xs, ys, (x, y) -> sol(x, y, t); camera=(15, 80), colorbar=false, color)
        plot!(; xguide="x", yguide="y", zguide="u")
        plot!(; zlim=(-0.01solmax, solmax))
        plot!(; size=(400, 400), margin=-20Plots.mm)
    end
    gif(anim, fn; fps)
end

function show_tau(u, v; disp=false)
    @syms x y t
    tau = make_soliton_tau(u, v)
    if disp
        print("tau(x, y, t) = ")
        display(tau(x, y, t))
    else
        @show tau(x, y, t)
    end
    
end

function show_tau_and_gif_sol(sol, u, v, xs, ys, ts;
        fn="tmp", fps=20, color=rescale_cgrad(:gist_earth),
        disp=false
    )
    show_tau(u, v; disp)
    flush(stdout)
    
    gif_sol_3d(sol, xs, ys, ts; fn, fps, color) |> display
    plot_sol_3d(sol, xs, ys, 0.0; color) |> display
    gif_sol_heatmap(sol, xs, ys, ts;
        fn=replace(fn, r"\.gif"=>"_heatmap.gif"), fps, color) |> display
    plot_sol_heatmap(sol, xs, ys, 0.0; color)
end

show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts;
        fn="tmp.gif", fps=20, color=rescale_cgrad(:gist_earth),
        disp=false) =
    show_tau_and_gif_sol(sol, (P-k)/2, (P+k)/2, xs, ys, ts; fn, fps, color, disp)

# %%
u, v = [1.0], [-0.8]

@syms x y t
tau = make_soliton_tau(u, v)
2SymPy.diff(log(tau(x, y, t)), x, x).factor() |> display
sol = make_soliton_solution(u, v)
sol(x, y, t).factor() |> display

sol = make_soliton_solution(u, v)
xs = range(-25, 25, 251)
ys = range(-25, 25, 251)
ts = range(-25, 25, 101)
show_tau_and_gif_sol(sol, u, v, xs, ys, ts; fn="1-soliton.gif", disp=true)

# %%
u, v = [1.0, 1.5], [-0.8, -1.0]
sol = make_soliton_solution(u, v)
xs = range(-25, 25, 251)
ys = range(-25, 25, 251)
ts = range(-25, 25, 101)
show_tau_and_gif_sol(sol, u, v, xs, ys, ts; fn="2-soliton.gif", disp=true)

# %%
k, P = [0.5, 0.5], [0.66, -0.66]
sol = make_soliton_solution_AB(k, P)
xs = range(-50, 50, 251)
ys = range(-50, 50, 251)
ts = range(-200, 200, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG1.gif", disp=true)

# %%
k, P = [0.5, 1.0], [0.75, 0.25]
sol = make_soliton_solution_AB(k, P)
xs = range(-50, 50, 251)
ys = range(-75, 75, 251)
ts = range(-200, 200, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG2.gif", disp=true)

# %%
k, P = [0.5, 0.5], [-0.26, 0.75]
sol = make_soliton_solution_AB(k, P)
xs = range(-50, 50, 251)
ys = range(-75, 75, 251)
ts = range(-200, 200, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG3.gif", disp=true)

# %%
k, P = [0.5, 0.5], [-0.5+1e-10, 0.5+1e-9]
sol = make_soliton_solution_AB(k, P)
xs = range(-50, 50, 251)
ys = range(-75, 75, 251)
ts = range(-200, 200, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG4.gif", disp=true)

# %%
k, P = [1.0, 0.5], [0.5-1e-7, 0.0]
sol = make_soliton_solution_AB(k, P)
xs = range(-100, 100, 251)
ys = range(-150, 150, 251)
ts = range(-400, 400, 101)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG5.gif", disp=true)

# %%
k, P = [1.0, 2.0, 3.0], [-0.333, -0.667, -1.667]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="AB_FIG6.gif")

# %%
k, P = [1.0, 2.0, 3.0], [-0.333, -0.667, -1.66]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_a.gif")

# %%
k, P = [1.0, 2.0, 3.0], [-0.333, -0.667, -1.5]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_b.gif")

# %%
k, P = [0.5, 1.0, 2.3], [0.75, 0.25, -1.0]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_c.gif")

# %%
k, P = [0.5, 1.0, 1.5], [0.75, 0.25, -0.25]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 50, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_d.gif")

# %%
k, P = [0.5, 1.0, 2.25], [0.75, 0.25, -1.0]
sol = make_soliton_solution_AB(k, P)
xs = range(-25, 25, 251)
ys = range(-25, 50, 251)
ts = range(-25, 25, 201)
show_tau_and_gif_sol_AB(sol, k, P, xs, ys, ts; fn="3-soliton_e.gif")

# %%
