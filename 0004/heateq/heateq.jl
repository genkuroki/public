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

n = 200
x = range(-π, π; length=n+1)[1:end-1]
dx = step(x)
m = 200
noise1, noise2 = 0.3randn(m), 0.3randn(m)
f(x) = SVector{m}((1 .+ noise1)*sin(x) .- (1 .+ noise2)*cos(2x))
u0 = f.(x);

using Plots, Printf
ENV["GKSwstype"]=100

P = plot(; legend=false)
for k in 1:m
    plot!(x, (p -> p[k]).(u0); lw=0.1, color=:blue)
end
plot()
savefig(P, "heateq0.png")

tmax = 1.0
t, u = heateq4(u0, dx, tmax);

ylim = (minimum(minimum.(u0)), maximum(maximum.(u0)))
anim = @animate for i in 1:length(t)
    title = @sprintf("t = %4.2f", t[i])
    plot(; legend=false)
    for k in 1:m
        plot!(x, (p -> p[k]).(u[:, i]); ylim, lw=0.1, color=:blue)
    end
    title!(title)
end
gif(anim, "heateq.gif")
