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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
@show VERSION
@time using DifferentialEquations
@time using Plots
@time using ProgressMeter
@time using Fmt # pkg> add https://github.com/bicycle1885/Fmt.jl

# %%
x ⪅ y = x < y || x ≈ y

function isin_ngon(n, R, x, y)
    α = 2π/n
    θ = mod(angle(x + im*y), α)
    r = √(x^2 + y^2) / R
    X, Y = r*cos(θ), r*sin(θ)
    (1 - cos(α))*Y ⪅ sin(α)*(1 - X)
end

isin_unit_pentagon(x, y) = isin_ngon(5, 1, x, y)

x = -1:0.1:1
y = -1:0.1:1
D = isin_unit_pentagon.(x, y')
heatmap(x, y, D'; colorbar=false, size=(200, 200))

# %%
module My

struct Param{TX, TY, TH, TV, TE, TB, TV_dic}
    x::TX
    y::TY
    h⁻¹::TH
    h⁻²::TH
    V::TV
    E::TE
    B::TB
    V_dic::TV_dic
end

Base.length(p::Param) = length(p.V)
Base.size(p::Param) = (length(p.x), length(p.y))

mystep(x::AbstractVector) = x[2] - x[1]
mystep(x::StepRangeLen) = step(x)

function Param(isin_func, x::AbstractVector, y::AbstractVector)
    @assert mystep(x) == mystep(y)
    h⁻¹ = 1/mystep(x)
    h⁻² = 1/mystep(x)^2
    I, J = eachindex(x), eachindex(y)
    D = isin_func.(x, y')
    V = Tuple{Int64, Int64}[]
    E_mat = [Tuple{Int64, Int64}[] for i in I, j in J]
    B_mat = [Tuple{Int64, Int64}[] for i in I, j in J]
    for j in J, i in I
        if D[i, j]
            push!(V, (i, j))
            for (p, q) in ((-1, 0), (1, 0), (0, -1), (0, 1))
                if i+p in I && j+q in J && D[i+p, j+q]
                    push!(E_mat[i, j], (i+p, j+q))
                else
                    D[i-p, j-q] && push!(B_mat[i, j], (i-p, j-q))
                end
            end
        end
    end
    V_dic = Dict(V[k] => k for k in eachindex(V))
    E = [getindex.(Ref(V_dic), E_mat[V[k]...]) for k in eachindex(V)]
    B = [getindex.(Ref(V_dic), B_mat[V[k]...]) for k in eachindex(V)]
    Param(x, y, h⁻¹, h⁻², V, E, B, V_dic)
end

function vec2mat!(p::Param, u, U)
    (; V) = p
    for k in eachindex(V)
        U[V[k]...] = u[k]
    end
    U
end
    
function vec2mat(p::Param, u; default_value=NaN)
    U = fill(eltype(u)(default_value), size(p))
    vec2mat!(p, u, U)
end

function mat2vec!(p::Param, U, u)
    (; V_dic) = p
    for t in keys(V_dic)
        u[V_dic[t]] = U[t...]
    end
    u
end

function mat2vec(p::Param, U)
    u = similar(eltype(U)[], length(p))
    mat2vec!(p, U, u)
end

end

x = -1:0.1:1
y = -1:0.1:1
p = My.Param(isin_unit_pentagon, x, y)

(; h⁻², V, E, B) = p
@show h⁻²

for k in (18, 19, 32, 131)
    @show k
    @show V[k]
    @show E[k]
    @show B[k]
    @show V[E[k]]
    @show V[B[k]]
    println("-"^80)
end

u0 = collect(1.0:length(V))
U0 = Int.(My.vec2mat(p, u0; default_value=0))
@show My.mat2vec(p, U0) == u0
U0'

# %%
function anim_sol2D(sol; L=201, fps=20, fn="tmp.gif",
        clim=(-0.06, 0.06), color=:bwr, bg=:lightgray)
    (; prob, u) = sol
    (; tspan, p) = prob
    (; x, y) = p
    ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
    prog = Progress(length(ts), 0)
    anim = @animate for t in ts
        U = My.vec2mat(p, sol(t).x[2])
        heatmap(x, y, U'; colorbar=false, clim, color, bg)
        title!(f"t = {$t:4.2f}"; titlefontsize=12)
        plot!(size=(420, 440), xlim=extrema(x), ylim=extrema(y), aspectratio=1)
        next!(prog)
    end
    gif(anim, fn; fps)
end

function plot_sol(sol, t; zlim=(-0.06, 0.06),
        color=:bwr, bg=:lightgray, camera=(30, 85), size=(640, 480), kwargs...)
    (; prob, u) = sol
    (; tspan, p) = prob
    (; x, y) = p
    U = My.vec2mat(p, sol(t).x[2])
    surface(x, y, U'; colorbar=false, zlim, color, bg, camera, size)
    title!(f"t = {$t:5.2f}"; titlefontsize=12)
    plot!(; kwargs...)
end

function anim_sol(sol; L=101, fps=10, fn="tmp.gif", zlim=(-0.06, 0.06),
        color=:bwr, bg=:lightgray, camera=(30, 85), size=(640, 480), kwargs...)
    (; prob, u) = sol
    (; tspan, p) = prob
    (; x, y) = p
    umin = minimum(minimum(u.x[2]) for u in u)
    umax = maximum(maximum(u.x[2]) for u in u)
    ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
    prog = Progress(length(ts), 0)
    anim = @animate for t in ts
        plot_sol(sol, t; color, bg, camera, size, zlim, kwargs...)
        next!(prog)
    end
    gif(anim, fn; fps)
end

# %%
function g!(dv, v, u, p, t)
    (; h⁻², E) = p
    @inbounds for k in eachindex(dv)
        dv[k] = h⁻² * sum(u[l] - u[k] for l in E[k])
    end
end

function g_dirichlet!(dv, v, u, p, t)
    (; h⁻², E) = p
    @inbounds for k in eachindex(dv)
        dv[k] = h⁻² * (
            sum(u[l] - u[k] for l in E[k]; init=0.0) - (4 - length(E[k])) * u[k]
        )
    end
end

function g_neumann!(dv, v, u, p, t)
    (; h⁻², E, B) = p
    @inbounds for k in eachindex(dv)
        dv[k] = h⁻² * (
            sum(u[l] - u[k] for l in E[k]; init=0.0) + sum(u[l] - u[k] for l in B[k]; init=0.0)
        )
    end
end

function g_absorbing!(dv, v, u, p, t)
    (; h⁻¹, h⁻², E, B) = p
    @inbounds for k in eachindex(dv)
        dv[k] = h⁻² * (
            sum(u[l] - u[k] for l in E[k]; init=0.0) + sum(u[l] - u[k] for l in B[k]; init=0.0)
        ) - h⁻¹ * 2length(B[k]) * v[k]
    end
end

# %%
x = y = range(-1, 1; length=21)
p = My.Param(isin_unit_pentagon, x, y)
tspan = (0.0, 2.0)

F(a, b, c, x, y) = exp(-((x-a)^2+(y-b)^2)/(2c^2))
f(x, y) = F(0.4, 0.1, 0.1, x, y)
U0 = f.(x, y')
V0 = zero(U0)

u0 = My.mat2vec(p, U0)
v0 = My.mat2vec(p, V0);

# %%
tspan = (0.0, 2.0)
prob = SecondOrderODEProblem(g_dirichlet!, v0, u0, tspan, p)

sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.01)

(; u) = sol
umin = minimum(minimum(u.x[2]) for u in u)
umax = maximum(maximum(u.x[2]) for u in u)
@show umin, umax;

# %%
anim_sol2D(sol; fn="wave eq on pentagon with dirichlet bc 2d Part 4.gif", clim=(-0.15, 0.15))

# %%
plot_sol(sol, 1.0; zlim=(-0.15, 0.15))

# %%
anim_sol(sol; fn="wave eq on pentagon with dirichlet bc Part 4.gif", zlim=(-0.15, 0.15))

# %%
tspan = (0.0, 2.0)
prob = SecondOrderODEProblem(g_neumann!, v0, u0, tspan, p)

sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.01)

(; u) = sol
umin = minimum(minimum(u.x[2]) for u in u)
umax = maximum(maximum(u.x[2]) for u in u)
@show umin, umax;

# %%
anim_sol2D(sol; fn="wave eq on pentagon with neumann 2d bc Part 4.gif", clim=(-0.15, 0.15))

# %%
plot_sol(sol, 1.0; zlim=(-0.15, 0.15))

# %%
anim_sol(sol; fn="wave eq on pentagon with neumann bc Part 4.gif", zlim=(-0.15, 0.15))

# %%
tspan = (0.0, 3.0)
prob = SecondOrderODEProblem(g_absorbing!, v0, u0, tspan, p)

sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.01)

(; u) = sol
umin = minimum(minimum(u.x[2]) for u in u)
umax = maximum(maximum(u.x[2]) for u in u)
@show umin, umax;

# %%
anim_sol2D(sol; L=301, fn="wave eq on pentagon with absorbing bc 2d Part 4.gif", clim=(-0.15, 0.15))

# %%
plot_sol(sol, 2.0; zlim=(-0.15, 0.15))

# %%
anim_sol(sol; L=151, fn="wave eq on pentagon with absorbing bc Part 4.gif", zlim=(-0.15, 0.15))

# %%
