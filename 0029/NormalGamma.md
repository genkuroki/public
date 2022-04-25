---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.3
  kernelspec:
    display_name: Julia 1.7.2
    language: julia
    name: julia-1.7
---

```julia
module O

using Distributions
using Random

struct NormalGamma <: ContinuousUnivariateDistribution
    μ
    λ
    a
    b
end

function Base.rand(rng::AbstractRNG, d::NormalGamma)
    (; μ, λ, a, b) = d
    rand(rng, Normal(μ, inv(sqrt(λ * rand(Gamma(a, b))))))
end

end
```

```julia
tmp = Vector{Float64}(undef, 10^6)
@time O.rand!(O.NormalGamma(1, 2, 2.1, 0.5), tmp)
@time O.rand!(O.NormalGamma(1, 2, 2.1, 0.5), tmp)
@time O.rand!(O.NormalGamma(1, 2, 2.1, 0.5), tmp);
```

```julia
@code_warntype rand(O.Random.default_rng(), O.NormalGamma(1, 2, 3, 4))
```

```julia
module Q

using Distributions
using Random

struct NormalGamma{Tμ, Tλ, Ta, Tb} <: ContinuousUnivariateDistribution
    μ::Tμ
    λ::Tλ
    a::Ta
    b::Tb
end

function Base.rand(rng::AbstractRNG, d::NormalGamma)
    (; μ, λ, a, b) = d
    rand(rng, Normal(μ, inv(sqrt(λ * rand(Gamma(a, b))))))
end

end
```

```julia
using Random
using Distributions
using StatsPlots
default(fmt = :png)

res(xlim, A) = A[first(xlim) .< A .< last(xlim)]

μ, λ, a, b = 1, 2, 1.5, 4
ng_org = O.NormalGamma(μ, λ, a, b)
ng_rev = Q.NormalGamma(μ, λ, a, b)
normal = Normal(μ, 1/√(λ*a*b))
tdist = μ + TDist(2a)/√(λ*a*b)

A = Vector{Float64}(undef, 10^6)
B = similar(A)
C = similar(A)
A = @time rand!(ng_org, A)
A = @time rand!(ng_org, A)
A = @time rand!(ng_org, A)
B = @time rand!(ng_rev, B)
B = @time rand!(ng_rev, B)
B = @time rand!(ng_rev, B)
C = @time rand!(tdist, C)
C = @time rand!(tdist, C)
C = @time rand!(tdist, C)
;
```

```julia
xlim = quantile.(Ref(A), (0.001, 0.999))
bin = round(Int, abs(-(xlim...))*50)
stephist(res(xlim, A); norm=true, label="original NormalGamma", bin)
stephist!(res(xlim, B); norm=true, label="revised NormalGamma", bin, ls=:dash)
stephist!(res(xlim, C); norm=true, label="μ+TDist(2a)/√(λ*a*b)", bin, ls=:dashdot)
plot!(normal, xlim...; label="normal approx.", c=:red, ls=:dot)
plot!(; xlim)
```

```julia

```
