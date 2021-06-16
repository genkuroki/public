using Pkg
pkgname = "HeatEqCalc"

if !isempty(ARGS)

@show ARGS

code_pkg = raw"""
module HeatEqCalc

function laplacian4!(v, u, dx)
    v[begin] = (u[end] + u[begin+1] - 2u[begin])/dx^2
    @. @views v[begin+1:end-1] = (u[begin:end-2] + u[begin+2:end] - 2u[begin+1:end-1])/dx^2
    v[end] = (u[end-1] + u[begin] - 2u[1end])/dx^2
    return
end

function heateq4(u0, dx, tmax, N=200)
    t = 0:dx:tmax
    dt = step(t)
    u = similar(u0, length(u0), length(t)+1)
    u[:, 1] = u0
    v = similar(u0)
    N = 100
    for i in 2:length(t)+1
        @. @views u[:, i] = u[:, i-1]
        for _ in 1:N
            @views laplacian4!(v, u[:, i], dx)
            @. @views u[:, i] += v*dt/N
        end
    end
    t, u
end

using StaticArrays

function calc_sol(
        n = 200,
        x = range(-π, π; length=n+1)[1:end-1],
        m = 200
    )
    dx = step(x)
    noise1, noise2 = 0.3randn(m), 0.3randn(m)
    f(x) = SVector{m}((1 .+ noise1)*sin(x) .- (1 .+ noise2)*cos(2x))
    u0 = f.(x)
    tmax = 1.0
    t, u = heateq4(u0, dx, tmax)
    (; t, x, u, m)
end

using DelimitedFiles

function save_sol(sol)
    t, x, u, m = sol
    U = vcat((getindex.(u, k) for k in 1:m)...)
    writedlm("sol_m.txt", m)
    writedlm("sol_t.txt", t)
    writedlm("sol_x.txt", x)
    writedlm("sol_U.txt", U)
end

function load_sol(name)
    m = parse(Int, read(name * "_m.txt", String))
    t = vec(readdlm(name * "_t.txt"))
    x = vec(readdlm(name * "_x.txt"))
    U = readdlm(name * "_U.txt")
    u = [SVector((U[k*length(x) + i, j] for k in 0:m-1)...) for i in 1:length(x), j in 1:length(t)]
    (; t, x, u, m)
end

if isfile(joinpath(pkgdir(@__MODULE__), "src", "trace.jl"))
    include("trace.jl")
end

end # module
"""

code_deps = raw"""

[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
"""

isdir(pkgname) || Pkg.generate(pkgname)

project_toml = joinpath(pkgname, "Project.toml")
code_toml = replace(read(project_toml, String), r"\[deps\](.|\n)*$"=>"")
write(project_toml, code_toml * code_deps)

pkg_jl = joinpath(pkgname, "src", pkgname * ".jl")
write(pkg_jl, code_pkg)

run(`julia --compile=all --trace-compile=$pkgname/src/trace.jl main.jl`)

Pkg.activate(pkgname)
using HeatEqCalc

else

Pkg.activate(pkgname)
using HeatEqCalc
sol = HeatEqCalc.calc_sol()
HeatEqCalc.save_sol(sol)

end