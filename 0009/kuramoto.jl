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
a = 1.3; tmax = 25.0; using Distributions, DifferentialEquations, Plots
kuramoto!(dθ, θ, param, t) = (N = length(θ); (K, ω) = param; 
    for i in 1:N dθ[i] = ω[i] + K*mean(sin(θ[j] - θ[i]) for j in 1:N) end)
(m, n) = (32, 16); (d, v) = (Normal(0, 2), 1.0); K_c = 2/π/pdf(d, 0)
θ₀ = 2π*rand(m, n); tspan = (0.0, tmax); K = a*K_c; ω = rand(d, m, n) .+ v
sol = solve(ODEProblem(kuramoto!, θ₀, tspan, (K, ω)))
anim = @animate for t in [fill(0, 20); 0:0.1:tmax; fill(tmax, 20)]
    plot(size=(400, 220), title="Kuramoto model K = $(a)K_c: t = $t", titlefontsize=10)
    heatmap!(sin.(sol(t))'; c=:bwr, colorbar=false, frame=false, ticks=false)
end; gif(anim, "kuramoto$(a).gif")

# %%
quote
a = 1.3; tmax = 25.0; using Distributions, DifferentialEquations, Plots
kuramoto!(dθ, θ, param, t) = (N = length(θ); (K, ω) = param; 
    for i in 1:N dθ[i] = ω[i] + K*mean(sin(θ[j] - θ[i]) for j in 1:N) end)
(m, n) = (32, 16); (d, v) = (Normal(0, 2), 1.0); K_c = 2/π/pdf(d, 0)
θ₀ = 2π*rand(m, n); tspan = (0.0, tmax); K = a*K_c; ω = rand(d, m, n) .+ v
sol = solve(ODEProblem(kuramoto!, θ₀, tspan, (K, ω)))
anim = @animate for t in [fill(0, 20); 0:0.1:tmax; fill(tmax, 20)]
    plot(size=(400, 220), title="Kuramoto model K = $(a)K_c: t = $t", titlefontsize=10)
    heatmap!(sin.(sol(t))'; c=:bwr, colorbar=false, frame=false, ticks=false)
end; gif(anim, "kuramoto$(a).gif")
end |> Base.remove_linenums! |> x -> display("text/markdown", "```julia\n$x\n```")

# %%
begin
    a = 1.3
    tmax = 25.0
    using Distributions, DifferentialEquations, Plots
    kuramoto!(dθ, θ, param, t) = begin
            N = length(θ)
            (K, ω) = param
            for i = 1:N
                dθ[i] = ω[i] + K * mean((sin(θ[j] - θ[i]) for j = 1:N))
            end
        end
    (m, n) = (32, 16)
    (d, v) = (Normal(0, 2), 1.0)
    K_c = (2 / π) / pdf(d, 0)
    θ₀ = (2π) * rand(m, n)
    tspan = (0.0, tmax)
    K = a * K_c
    ω = rand(d, m, n) .+ v
    sol = solve(ODEProblem(kuramoto!, θ₀, tspan, (K, ω)))
    anim = @animate(for t = [fill(0, 20); 0:0.1:tmax; fill(tmax, 20)]
                plot(size = (400, 220), title = "Kuramoto model K = $(a)K_c: t = $(t)", titlefontsize = 10)
                heatmap!((sin.(sol(t)))'; c = :bwr, colorbar = false, frame = false, ticks = false)
            end)
    gif(anim, "kuramoto$(a).gif")
end

# %%
