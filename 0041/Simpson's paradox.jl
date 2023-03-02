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
#     display_name: Julia 1.9.0-beta4
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6, guidefontsize=10, legendfontsize=10)

# %%
function randxyz(; distz, distx, disty)
    z = rand(distz)
    x = rand(distx(z))
    y = rand(disty(x, z))
    x, y, z
end

function randxyz(n ;
        distz = 5 + Normal(), 
        distx = z -> z + Normal(),
        disty = (x, z) -> -x + 4z + Normal(),
    )
    stack(collect(randxyz(; distz, distx, disty)) for _ in 1:n)
end

# %%
Random.seed!(4649373)
color = :CMRmap

n = 1000
xyz = randxyz(n)
x, y, z = eachrow(xyz)

@show fit(MvNormal, xyz)

@show fit(Normal, z)
@show c = [ones(n) z] \ x
@show fit(Normal, x - (c[1] .+ c[2]*z))
@show [ones(n) x z] \ y
println()
@show fit(Normal, x)
@show c = [ones(n) x] \ z
@show fit(Normal, z - (c[1] .+ c[2]*x))

scatter(x, y; label="(xᵢ, yᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
plot!(xguide="xᵢ", yguide="yᵢ", colorbar_title="round(zᵢ)")

# %%
A = [ones(n) x]
@show a = A \ y

scatter(x, y; label="(xᵢ, yᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
plot!(xguide="xᵢ", yguide="yᵢ", colorbar_title="zᵢ")
plot!(x -> a[1] + a[2]*x; label="", c=:red)

# %%
ŷ = @. a[1] + a[2]*x

scatter(x, y - ŷ; label="(xᵢ, yᵢ - ŷᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
plot!(xguide="xᵢ", yguide="yᵢ - ŷᵢ", colorbar_title="zᵢ")
hline!([0]; label="", c=:red)

# %%
ŷ = @. a[1] + a[2]*x

perm = sortperm(x)
scatter(1:n, (y - ŷ)[perm]; label="", msw=0, alpha=0.5, ms=3, marker_z=round.(z[perm]), color)
plot!(xguide="rank of xᵢ", yguide="yᵢ - ŷᵢ", colorbar_title="zᵢ")
hline!([0]; label="", c=:red)

# %%
B = [ones(n) x z]
@show b = B \ y

scatter(x, y; label="(xᵢ, yᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
for z in range(round.(extrema(z))...)
    plot!(x -> b[1] + b[2]*x + b[3]*z, z-2, z+2; label="", c=:blue)
end
plot!()

# %%
scatter(x, y - b[2]*x; label="(xᵢ, yᵢ - b₂xᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
for z in range(round.(extrema(z))...)
    plot!(x -> b[1] + b[3]*z, z-2, z+2; label="", c=:blue)
end
plot!()

# %%
function randxyz2(; distx, distz, disty)
    x = rand(distx)
    z = rand(distz(x))
    y = rand(disty(x, z))
    x, y, z
end

function randxyz2(n ;
        distx = 5 + Normal(0, √2), 
        distz = x -> 2.5 + 0.5x + Normal(0, 1/√2),
        disty = (x, z) -> -x + 4z + Normal(),
    )
    stack(collect(randxyz2(; distx, distz, disty)) for _ in 1:n)
end

# %%
n = 1000
xyz = randxyz2(n)
x, y, z = eachrow(xyz)

@show fit(MvNormal, xyz)

@show fit(Normal, z)
@show c = [ones(n) z] \ x
@show fit(Normal, x - (c[1] .+ c[2]*z))
@show [ones(n) x z] \ y
println()
@show fit(Normal, x)
@show c = [ones(n) x] \ z
@show fit(Normal, z - (c[1] .+ c[2]*x))

scatter(x, y; label="(xᵢ, yᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
plot!(xguide="xᵢ", yguide="yᵢ", colorbar_title="round(zᵢ)")

# %%
A = [ones(n) x]
@show a = A \ y

scatter(x, y; label="(xᵢ, yᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
plot!(xguide="xᵢ", yguide="yᵢ", colorbar_title="zᵢ")
plot!(x -> a[1] + a[2]*x; label="", c=:red)

# %%
ŷ = @. a[1] + a[2]*x

scatter(x, y - ŷ; label="(xᵢ, yᵢ - ŷᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
plot!(xguide="xᵢ", yguide="yᵢ - ŷᵢ", colorbar_title="zᵢ")
hline!([0]; label="", c=:red)

# %%
ŷ = @. a[1] + a[2]*x

perm = sortperm(x)
scatter(1:n, (y - ŷ)[perm]; label="", msw=0, alpha=0.5, ms=3, marker_z=round.(z[perm]), color)
plot!(xguide="rank of xᵢ", yguide="yᵢ - ŷᵢ", colorbar_title="zᵢ")
hline!([0]; label="", c=:red)

# %%
B = [ones(n) x z]
@show b = B \ y

scatter(x, y; label="(xᵢ, yᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
for z in range(round.(extrema(z))...)
    plot!(x -> b[1] + b[2]*x + b[3]*z, z-2, z+2; label="", c=:blue)
end
plot!()

# %%
scatter(x, y - b[2]*x; label="(xᵢ, yᵢ - b₂xᵢ)", msw=0, alpha=0.5, ms=3, marker_z=round.(z), color)
for z in range(round.(extrema(z))...)
    plot!(x -> b[1] + b[3]*z, z-2, z+2; label="", c=:blue)
end
plot!()

# %%
