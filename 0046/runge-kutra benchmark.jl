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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# * https://twitter.com/kame_no_mori/status/1734888762914820124
# * https://twitter.com/genkuroki/status/1739876470682697877
#
# >\#Julia言語 https://kamemori.com/research/fortran/speed_harmonic_oscillator_1d_ja.html のJuliaのコードに18文字追加して比較し直してみました。3回実行して最短の値を採用。<br>
# ><br>
# >1.64秒 ← gfortran -O3<br>
# >3.50秒 ← julia newton.jl<br>
# >1.63秒 ← julia newton_inline_fastmath.jl (18文字追加)<br>
# >1.18秒 ← julia> include("newton_inline_fastmath.jl")<br>
#
# * https://x.com/genkuroki/status/1739999710998618530?s=20
#
# >#Julia言語 さらなる訂正<br>
# ><br>
# >Julia版では34文字追加しないとかなり遅くなるという結果<br>
# ><br>
# >1.64秒 ← gfortran -O3<br>
# >3.50秒 ← julia newton.jl<br>
# >1.41秒 ← julia newton_inline_fastmath_inbounds_simd.jl (34文字追加)<br>
# >1.06秒 ← julia> include("newton_inline_fastmath.jl")<br>

# %% [markdown]
# __gfortranとの比較__

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/runge-kutta%20benchmark/2023-12-27a%20gfortran%20-O3.png">

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/runge-kutta%20benchmark/2023-12-27b%20julia%20newton.jl.png">

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/runge-kutta%20benchmark/2023-12-27c%20julia%20newton_inline_fastmath.jl.png">

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/runge-kutta%20benchmark/2023-12-27d%20julia%20include.png">

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/runge-kutta%20benchmark/2023-12-27e%20julia%20newton_inline_fastmath_inbounds_simd.jl.png">

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/runge-kutta%20benchmark/2023-12-27f%20julia%20include%20newton_inline_fastmath_inbounds_simd.jl.png">

# %%
function main(nt = 100000000)
    mass = 1.0
    k = 1.0
    dt = 1e-2

    xt = zeros(Float64, nt+1)
    vt = zeros(Float64, nt+1)

    x = 0.0
    v = 1.0

    for it = 1:nt+1
        xt[it] = x
        vt[it] = v
        x, v = Runge_Kutta_4th(x, v, dt, mass, k)
    end

    open("result_julia.out", "w") do file
        for it = nt-999:nt
            println(file, "$(it*dt) $(xt[it]) $(vt[it])")
        end
    end
end

function Runge_Kutta_4th(x, v, dt, mass, k)
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

@time main()
@time main()
@time main()

# %%
function main_inline(nt = 100000000)
    mass = 1.0
    k = 1.0
    dt = 1e-2

    xt = zeros(Float64, nt+1)
    vt = zeros(Float64, nt+1)

    x = 0.0
    v = 1.0

    for it = 1:nt+1
        xt[it] = x
        vt[it] = v
        x, v = Runge_Kutta_4th_inline(x, v, dt, mass, k)
    end

    open("result_julia.out", "w") do file
        for it = nt-999:nt
            println(file, "$(it*dt) $(xt[it]) $(vt[it])")
        end
    end
end

@inline function Runge_Kutta_4th_inline(x, v, dt, mass, k)
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

@time main_inline()
@time main_inline()
@time main_inline()

# %%
function main_inline_simd(nt = 100000000)
    mass = 1.0
    k = 1.0
    dt = 1e-2

    xt = zeros(Float64, nt+1)
    vt = zeros(Float64, nt+1)

    x = 0.0
    v = 1.0

    @simd for it in 1:nt+1
        xt[it] = x
        vt[it] = v
        x, v = Runge_Kutta_4th_inline(x, v, dt, mass, k)
    end

    open("result_julia.out", "w") do file
        for it = nt-999:nt
            println(file, "$(it*dt) $(xt[it]) $(vt[it])")
        end
    end
end

@inline function Runge_Kutta_4th_inline(x, v, dt, mass, k)
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

@time main_inline_simd()
@time main_inline_simd()
@time main_inline_simd()

# %%
function main_fastmath(nt = 100000000)
    mass = 1.0
    k = 1.0
    dt = 1e-2

    xt = zeros(Float64, nt+1)
    vt = zeros(Float64, nt+1)

    x = 0.0
    v = 1.0

    for it in 1:nt+1
        xt[it] = x
        vt[it] = v
        x, v = Runge_Kutta_4th_fastmath(x, v, dt, mass, k)
    end

    open("result_julia.out", "w") do file
        for it = nt-999:nt
            println(file, "$(it*dt) $(xt[it]) $(vt[it])")
        end
    end
end

@fastmath function Runge_Kutta_4th_fastmath(x, v, dt, mass, k)
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

@time main_fastmath()
@time main_fastmath()
@time main_fastmath()

# %%
function main_inline_fastmath(nt = 100000000)
    mass = 1.0
    k = 1.0
    dt = 1e-2

    xt = zeros(Float64, nt+1)
    vt = zeros(Float64, nt+1)

    x = 0.0
    v = 1.0

    for it in 1:nt+1
        xt[it] = x
        vt[it] = v
        x, v = Runge_Kutta_4th_inline_fastmath(x, v, dt, mass, k)
    end

    open("result_julia.out", "w") do file
        for it = nt-999:nt
            println(file, "$(it*dt) $(xt[it]) $(vt[it])")
        end
    end
end

@inline @fastmath function Runge_Kutta_4th_inline_fastmath(x, v, dt, mass, k)
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

@time main_inline_fastmath()
@time main_inline_fastmath()
@time main_inline_fastmath()

# %%
function main_inline_fastmath_inbounds(nt = 100000000)
    mass = 1.0
    k = 1.0
    dt = 1e-2

    xt = zeros(Float64, nt+1)
    vt = zeros(Float64, nt+1)

    x = 0.0
    v = 1.0

    @inbounds for it in 1:nt+1
        xt[it] = x
        vt[it] = v
        x, v = Runge_Kutta_4th_inline_fastmath(x, v, dt, mass, k)
    end

    open("result_julia.out", "w") do file
        for it = nt-999:nt
            println(file, "$(it*dt) $(xt[it]) $(vt[it])")
        end
    end
end

@inline @fastmath function Runge_Kutta_4th_inline_fastmath(x, v, dt, mass, k)
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

@time main_inline_fastmath_inbounds()
@time main_inline_fastmath_inbounds()
@time main_inline_fastmath_inbounds()

# %%
function main_inline_fastmath_inbounds_simd(nt = 100000000)
    mass = 1.0
    k = 1.0
    dt = 1e-2

    xt = zeros(Float64, nt+1)
    vt = zeros(Float64, nt+1)

    x = 0.0
    v = 1.0

    @inbounds @simd for it in 1:nt+1
        xt[it] = x
        vt[it] = v
        x, v = Runge_Kutta_4th_inline_fastmath(x, v, dt, mass, k)
    end

    open("result_julia.out", "w") do file
        for it = nt-999:nt
            println(file, "$(it*dt) $(xt[it]) $(vt[it])")
        end
    end
end

@inline @fastmath function Runge_Kutta_4th_inline_fastmath(x, v, dt, mass, k)
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

@time main_inline_fastmath_inbounds_simd()
@time main_inline_fastmath_inbounds_simd()
@time main_inline_fastmath_inbounds_simd()

# %%
using BenchmarkTools

nt = 10^6
println("nt = ", nt)
print("main(nt):                              "); @btime main(nt)
print("main_inline(nt):                       "); @btime main_inline(nt)
print("main_inline_simd(nt):                  "); @btime main_inline_simd(nt)
print("main_fastmath(nt):                     "); @btime main_fastmath(nt)
print("main_inline_fastmath(nt):              "); @btime main_inline_fastmath(nt)
print("main_inline_fastmath_inbounds(nt):     "); @btime main_inline_fastmath_inbounds(nt)
print("main_inline_fastmath_inbounds_simd(nt):"); @btime main_inline_fastmath_inbounds_simd(nt)

# %%
versioninfo()

# %%
