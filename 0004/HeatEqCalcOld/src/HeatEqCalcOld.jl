module HeatEqCalcOld

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

save_sol(sol) = write("sol.txt", repr(sol))

end # module
