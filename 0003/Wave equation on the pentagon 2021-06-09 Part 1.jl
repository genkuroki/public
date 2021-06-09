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

# %%
@show VERSION
@time using DifferentialEquations
@time using Plots
@time using ProgressMeter
@time using Parameters

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

using Parameters: @unpack

struct Param{TX, TY, TH, TV, TE, TV_dic}
    x::TX
    y::TY
    h⁻²::TH
    V::TV
    E::TE
    V_dic::TV_dic
end

Base.length(p::Param) = length(p.V)
Base.size(p::Param) = (length(p.x), length(p.y))

mystep(x::AbstractVector) = x[2] - x[1]
mystep(x::StepRangeLen) = step(x)

function Param(isin_func, x::AbstractVector, y::AbstractVector)
    @assert mystep(x) == mystep(y)
    h⁻² = 1/mystep(x)^2
    I, J = eachindex(x), eachindex(y)
    D = isin_func.(x, y')
    V = Tuple{Int64, Int64}[]
    E_mat = [Tuple{Int64, Int64}[] for i in I, j in J]
    B_mat = zeros(Int64, length(I), length(J))
    for j in J, i in I
        if D[i, j]
            push!(V, (i, j))
            for q in (-1, 1), p in (-1, 1)
                if i+p in I && j+q in J && D[i+p, j+q]
                    push!(E_mat[i, j], (i+p, j+q))
                end
            end
        end
    end
    V_dic = Dict(V[k] => k for k in eachindex(V))
    E = [getindex.(Ref(V_dic), E_mat[V[k]...]) for k in eachindex(V)]
    Param(x, y, h⁻², V, E, V_dic)
end

function vec2mat!(p::Param, u, U)
    @unpack V = p
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
    @unpack V_dic = p
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

@unpack h⁻², V, E = p
@show h⁻²
@show V[50]
@show E[50]
@show V[E[50]]

u0 = collect(1.0:length(V))
U0 = My.vec2mat(p, u0; default_value=0)
@show My.mat2vec(p, U0) == u0
U0'

# %%
function anim_sol(sol; L=201, fps=20, fn="tmp.gif", color=:CMRmap, camera=(30, 60))
    @unpack prob, u = sol
    @unpack tspan, p = prob
    @unpack x, y = p
    umin = minimum(minimum(u.x[2]) for u in u)
    umax = maximum(maximum(u.x[2]) for u in u)
    ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
    prog = Progress(length(ts), 0)
    anim = @animate for t in ts
        U = My.vec2mat(p, sol(t).x[2])
        surface(x, y, U'; zlim=(umin, umax), colorbar=false, color, camera)
        title!("t = $(round(t, digits=1))"; titlefontsize=12)
        next!(prog)
    end
    gif(anim, fn; fps)
end

# %%
#f!(du, v, u, p, t) = @. du = v

function g_neumann!(dv, v, u, p, t)
    @unpack h⁻², E = p
    @inbounds for k in eachindex(dv)
        dv[k] = h⁻² * sum(u[l] - u[k] for l in E[k])
    end
end

x = y = range(-1, 1; length=401)
p = My.Param(isin_unit_pentagon, x, y)
tspan = (0.0, 20.0)

f(x, y) = exp(-(x^2 + y^2)/2)
U0 = f.(x, y')
V0 = zero(U0)

u0 = My.mat2vec(p, U0)
v0 = My.mat2vec(p, V0)

#@time prob = DynamicalODEProblem(g_neumann!, f!, v0, u0, tspan, p)
@time prob = SecondOrderODEProblem(g_neumann!, v0, u0, tspan, p)
sol_n = nothing; GC.gc()
@time sol_n = solve(prob; saveat=0.1)
sol_n = nothing; GC.gc()
@time sol_n = solve(prob; saveat=0.1);

# %%
anim_sol(sol_n; fps=10, fn="wave eq on pentagon neumann.gif")

# %%
#f!(du, v, u, p, t) = @. du = v

function g_dirichlet!(dv, v, u, p, t)
    @unpack h⁻², E = p
    @inbounds for k in eachindex(dv)
        dv[k] = h⁻² * (sum(u[l] - u[k] for l in E[k]) - (4 - length(E[k])) * u[k])
    end
end

x = y = range(-1, 1; length=401)
p = My.Param(isin_unit_pentagon, x, y)
tspan = (0.0, 20.0)

f(x, y) = exp(-(x^2 + y^2)/(2*0.2^2))
U0 = f.(x, y')
V0 = zero(U0)

u0 = My.mat2vec(p, U0)
v0 = My.mat2vec(p, V0)

#@time prob = DynamicalODEProblem(g_dirichlet!, f!, v0, u0, tspan, p)
@time prob = SecondOrderODEProblem(g_dirichlet!, v0, u0, tspan, p)
@time sol_d = solve(prob; saveat=0.1);

# %%
anim_sol(sol_d; fps=10, fn="wave eq on pentagon dirichlet CMRmap.gif")

# %%
anim_sol(sol_d; fps=10, fn="wave eq on pentagon dirichlet gist_earch.gif", color=:gist_earth)

# %%
