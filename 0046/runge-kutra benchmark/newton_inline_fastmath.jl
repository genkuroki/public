function main()
    mass = 1.0
    k = 1.0
    dt = 1e-2
    nt = 100000000

    xt = zeros(Float64, nt+1)
    vt = zeros(Float64, nt+1)

    x = 0.0
    v = 1.0

    for it = 1:nt+1
        xt[it] = x
        vt[it] = v
        x, v = Runge_Kutta_4th!(x, v, dt, mass, k)
    end

    open("result_julia.out", "w") do file
        for it = nt-999:nt
            println(file, "$(it*dt) $(xt[it]) $(vt[it])")
        end
    end
end

@inline @fastmath function Runge_Kutta_4th!(x, v, dt, mass, k)
    x1 = v
    v1 = force(x, mass, k)

    x2 = v + 0.5 * dt * v1
    v2 = force(x + 0.5 * x1 * dt, mass, k)

    x3 = v + 0.5 * dt * v2
    v3 = force(x + 0.5 * x2 * dt, mass, k)

    x4 = v + dt * v3
    v4 = force(x + x3 * dt, mass, k)

    x += (x1 + 2 * x2 + 2 * x3 + x4) * dt / 6
    v += (v1 + 2 * v2 + 2 * v3 + v4) * dt / 6

		return x, v
end

function force(x, mass, k)
    return -x * k / mass
end

main()
