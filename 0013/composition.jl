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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
module Q

using InteractiveUtils: subtypes

abstract type Creature end

Base.@kwdef mutable struct _Common
    hp::Int = 200
    atk::Int = 100
end

Base.@kwdef mutable struct Goblin <: Creature
    _common::_Common = _Common(88, 77)
    weapon::String = "dagger"
end

Base.@kwdef mutable struct Slime <: Creature
    _common::_Common = _Common(55, 44)
    sp::String = "acid"
end

for T in subtypes(Creature)
    n = length(fieldnames(_Common))
    @eval $(nameof(T))(a...) = $T(_Common(a[1:$n]...), a[$n+1:end]...)
    @eval Base.getproperty(x::$T, p::Symbol) =
        p ∈ fieldnames(_Common) ? getfield(x._common, p) : getfield(x, p)
    @eval Base.setproperty!(x::$T, p::Symbol, v) =
        p ∈ fieldnames(_Common) ? setfield!(x._common, p, v) : setfield!(x, p, v)
    @eval Base.propertynames(x::$T) = (fieldnames(_Common)...,
        (p for p in fieldnames($T) if first(string(p)) != '_')...)
    @eval function Base.show(io::IO, x::$T)
        props = getproperty.(Ref(x), propertynames(x))
        print(io, string(nameof($T)), '(', repr(first(props)))
        for p in props[2:end] print(io, ", ", repr(p)) end
        print(io, ')')
    end
end

end

@show Q.Goblin()
@show propertynames(Q.Goblin())
@show Q.Slime()
@show propertynames(Q.Slime())
println()
@show g = Q.Goblin(123, 45, "dirty dagger")
@show g.hp
@show g.atk
@show g.weapon
@show g.hp = 99
@show g
@show g.atk = 11
@show g
@show g.weapon = "nothing"
@show g
println()
@show s = Q.Slime(123, 45, "weak acid")
@show s.hp
@show s.atk
@show s.sp
@show s.hp = 66
@show s
@show s.atk = 99
@show s
@show s.sp = "strong acid"
@show s;

# %%
T = Q.Goblin
n = 2
:($(nameof(T))(a...) = $T(_Common(a[1:$n]...), a[$n+1:end]...)) |> Base.remove_linenums! |> print

# %%
:(Base.getproperty(x::$T, p::Symbol) =
    p ∈ fieldnames(_Common) ? getfield(x._common, p) : getfield(x, p)) |> Base.remove_linenums! |> print

# %%
:(Base.setproperty!(x::$T, p::Symbol, v) =
    p ∈ fieldnames(_Common) ? setfield!(x._common, p, v) : setfield!(x, p, v)) |> Base.remove_linenums! |> print

# %%
:(Base.propertynames(x::$T) = (fieldnames(_Common)...,
    (p for p in fieldnames($T) if first(string(p)) != '_')...)) |> Base.remove_linenums! |> print

# %%
:(function Base.show(io::IO, x::$T)
    props = getproperty.(Ref(x), propertynames(x))
    print(io, string(nameof($T)), '(', repr(first(props)))
    for p in props[2:end] print(io, ", ", repr(p)) end
    print(io, ')')
end) |> Base.remove_linenums! |> print

# %%
module O

using InteractiveUtils: subtypes

abstract type Creature end

Base.@kwdef mutable struct _Common
    hp::Int = 200
    atk::Int = 100
end
for p in fieldnames(_Common)
    @eval $p(x::Creature) = getfield(x._common, Symbol($p))
    @eval $(Symbol(:set_, p, :!))(x::Creature, v) = setfield!(x._common, Symbol($p), v)
end

Base.@kwdef mutable struct Goblin <: Creature
    _common::_Common = _Common(88, 77)
    weapon::String = "dagger"
end

Base.@kwdef mutable struct Slime <: Creature
    _common::_Common = _Common(55, 44)
    sp::String = "acid"
end

for T in subtypes(Creature)
    n = length(fieldnames(_Common))
    @eval $(nameof(T))(a...) = $T(_Common(a[1:$n]...), a[$n+1:end]...)
    for p in fieldnames(T)
        first(string(p)) == '_' && continue
        @eval $p(x::$T) = getfield(x, Symbol($p))
        @eval $(Symbol(:set_, p, :!))(x::$T, v) = setfield!(x, Symbol($p), v)
    end
end

end

@show O.Goblin()
@show O.Slime()
println()
@show g = O.Goblin(123, 45, "dirty dagger")
@show O.hp(g)
@show O.atk(g)
@show O.weapon(g)
@show O.set_hp!(g, 99)
@show g
@show O.set_atk!(g, 11)
@show g
@show O.set_weapon!(g, "nothing")
@show g
println()
@show s = O.Slime(123, 45, "weak acid")
@show O.hp(s)
@show O.atk(s)
@show O.sp(s)
@show O.set_hp!(s, 66)
@show s
@show O.set_atk!(s, 99)
@show s
@show O.set_sp!(s, "strong acid")
@show s;

# %%
