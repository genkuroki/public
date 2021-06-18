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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://discourse.julialang.org/t/ode-solvers-why-is-matlab-ode45-uncannily-stable/63052/15

# %%
using DifferentialEquations
using Plots

# %%
function wk5!(dP, P, params, t)

    #=

    The 5-Element WK with serial L
    set Ls = 0 for 4-Element WK parallel

    Formulation for DifferentialEquations.jl

    P:          solution vector (pressures p1 and p2)
    params: parameter tuple
                (Rc, Rp, C, Lp, Ls, I, q)


    I need to find a way to tranfer the function name as well
    for the time being we have to have the function in "I"

    =#

    # Split parameter tuple:
    Rc, Rp, C, Lp, Ls, I, q = params

    dP[1] = (
        -Rc / Lp * P[1]
        + (Rc / Lp - 1 / Rp / C) * P[2]
        + Rc * (1 + Ls / Lp) * didt(I, t, q)
        + I(t, q) / C
        )

    dP[2] = -1 / Rp / C * P[2] + I(t, q) / C

    return

end

# %%
# Generic Input Waveform
# max volume flow in ml/s
max_i = 425

# min volume flow in m^3/s
min_i = 0.0

T = 0.9

# Syst. Time in s
systTime = 2 / 5 * T

# Dicrotic notch time
dicrTime = 0.02

q_generic = (max_i, min_i, T, systTime, dicrTime)

function I_generic(t, q_generic)
    max_i, min_i, T, systTime, dicrTime = q_generic
    # implicit conditional using boolean multiplicator
    # sine waveform
    (
        (max_i - min_i) * sin(pi / systTime * (t % T))
        * (t % T < (systTime + dicrTime) )
        + min_i
    )
end

# %%
function didt(I, t, q)
    dt = 1e-3;
    didt = (I(t+dt, q) - I(t-dt, q)) / (2 * dt)
    return didt
end

# %%
plot(range(0, 2; length=2000), t -> didt(I_generic, t, q_generic); label="didt")

# %%
# Initial condition and time span
P0 = [0.0, 0.0]
tspan = (0.0, 30.0)

# Set parameters for Windkessel Model
Rc = 0.033
Rp = 0.6
C = 1.25
# L for serial model!
Ls = 0.01
# L for parallel
Lp = 0.02

I = I_generic
q = q_generic

p5 = (Rc, Rp, C, Lp, Ls, I, q)

problem = ODEProblem(wk5!, P0, tspan, p5)

dtmax = 1e-4

@time solutionTsit = solve(problem); GC.gc()
@time solutionTsitLowdt = solve(problem, Tsit5(), dtmax=dtmax); GC.gc()
@time solutionBS3 = solve(problem, BS3()); GC.gc()
@time solutionBS3Lowdt = solve(problem, BS3(), dtmax=dtmax); GC.gc()
@time solutionDP5 = solve(problem,DP5()); GC.gc()
@time solutionDP5Lowdt = solve(problem,DP5(), dtmax=dtmax); GC.gc()
@time solutionStiff = solve(problem, alg_hints=[:stiff]); GC.gc()

# %%
@time solutionTsit = solve(problem); GC.gc()
@time solutionTsitLowdt = solve(problem, Tsit5(), dtmax=dtmax); GC.gc()
@time solutionBS3 = solve(problem, BS3()); GC.gc()
@time solutionBS3Lowdt = solve(problem, BS3(), dtmax=dtmax); GC.gc()
@time solutionDP5 = solve(problem,DP5()); GC.gc()
@time solutionDP5Lowdt = solve(problem,DP5(), dtmax=dtmax); GC.gc()
@time solutionStiff = solve(problem, alg_hints=[:stiff]); GC.gc()

# %%
@show length(solutionTsit.t)
@show length(solutionTsitLowdt.t)
@show length(solutionBS3.t)
@show length(solutionBS3Lowdt.t)
@show length(solutionDP5.t)
@show length(solutionDP5Lowdt.t)
@show length(solutionStiff.t);

# %%
a, b = 0, 2

plot()
#plot!(t -> solutionTsit(t; idxs=1), a, b; label="Tsit")
plot!(t -> solutionTsitLowdt(t; idxs=1), a, b; label="TsitLowdt")
#plot!(t -> solutionBS3(t; idxs=1), a, b; label="BS3", ls=:dash)
plot!(t -> solutionBS3Lowdt(t; idxs=1), a, b; label="BS3Lowdt", ls=:dash)
#plot!(t -> solutionDP5(t; idxs=1), a, b; label="DP5", ls=:dashdot)
plot!(t -> solutionDP5Lowdt(t; idxs=1), a, b; label="DP5Lowdt", ls=:dashdot)
plot!(t -> solutionStiff(t; idxs=1), a, b; label="Stiff", ls=:dot, lw=1.5)

# %%
a, b = 0, 2

plot()
plot!(t -> solutionTsit(t; idxs=1), a, b; label="Tsit")
#plot!(t -> solutionTsitLowdt(t; idxs=1), a, b; label="TsitLowdt")
plot!(t -> solutionBS3(t; idxs=1), a, b; label="BS3", ls=:dash)
#plot!(t -> solutionBS3Lowdt(t; idxs=1), a, b; label="BS3Lowdt", ls=:dash)
plot!(t -> solutionDP5(t; idxs=1), a, b; label="DP5", ls=:dashdot)
#plot!(t -> solutionDP5Lowdt(t; idxs=1), a, b; label="DP5Lowdt", ls=:dashdot)
plot!(t -> solutionStiff(t; idxs=1), a, b; label="Stiff", ls=:dot, lw=1.5)

# %%
a, b = 28, 30

plot()
plot!(t -> solutionTsit(t; idxs=1), a, b; label="Tsit")
#plot!(t -> solutionTsitLowdt(t; idxs=1), a, b; label="TsitLowdt")
plot!(t -> solutionBS3(t; idxs=1), a, b; label="BS3", ls=:dash)
#plot!(t -> solutionBS3Lowdt(t; idxs=1), a, b; label="BS3Lowdt", ls=:dash)
plot!(t -> solutionDP5(t; idxs=1), a, b; label="DP5", ls=:dashdot)
#plot!(t -> solutionDP5Lowdt(t; idxs=1), a, b; label="DP5Lowdt", ls=:dashdot)
plot!(t -> solutionStiff(t; idxs=1), a, b; label="Stiff", ls=:dot, lw=1.5)

# %%
