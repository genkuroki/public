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
#     display_name: Julia 1.7.0-beta3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Statistics
using Plots

function mcpi1(n)
    c = 0
    for _ in 1:n
        c += rand()^2 + rand()^2 ≤ 1
    end
    4c/n
end

function mcpi2(n)
    c = 0
    for _ in 1:n
        x, y = rand(), rand()
        c += (x^2 + y^2 ≤ 1 && (1 - x)^2 + (1 - y)^2 ≤ 1)
    end
    2c/n + 2
end

function mcpi3(n)
    c = 0
    for _ in 1:n
        x, y = rand(), (1 - √3/2)*rand() + √3/2
        c += (x^2 + y^2 > 1 && (1 - x)^2 + y^2 > 1)
    end
    6(1 - √3/4 - (1 - √3/2)*c/n)
end

function sim(mcpi, n, L)
    X = Vector{Float64}(undef, L)
    Threads.@threads for i in 1:L
        X[i] = mcpi(n)
    end
    X
end

@show Threads.nthreads()
@time pi1 = sim(mcpi1, 10^6, 10^4)
@time pi2 = sim(mcpi2, 10^6, 10^4)
@time pi3 = sim(mcpi3, 10^6, 10^4)

@show mean(pi1), std(pi1)
@show mean(pi2), std(pi2)
@show mean(pi3), std(pi3)

plot()
stephist!(pi1; norm=true, label="mcpi1", lw=1.5, ls=:solid)
stephist!(pi2; norm=true, label="mcpi2", lw=1.5, ls=:dash)
stephist!(pi3; norm=true, label="mcpi3", lw=1.5, ls=:dashdot)

# %%
using Base64
showimg(mime, fn; tag="img") = open(fn) do f
    base64 = base64encode(f)
    display("text/html", """<$tag src="data:$mime;base64,$base64" />""")
end
showimg("image/jpeg", "1BD08057-D77D-44FF-9D3A-7620DA4CF740.jpeg"; tag="img width=80%")

# %%
