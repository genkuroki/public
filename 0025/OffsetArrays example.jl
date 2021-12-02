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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# https://twitter.com/dolphin7473/status/1466376598853148680

# %%
module O

using OffsetArrays

mutable struct Variables
    u::OffsetArray{ComplexF64, 1}
end

function update_u!(var::Variables)
    rhs = OffsetArray{ComplexF64}(undef, -2:2)
    rhs[-2:2] = randn(ComplexF64, 5)
    var.u .+= rhs
end

function main()
    u = OffsetArray{ComplexF64}(undef, -2:2)   
    var = Variables(u)
    var.u[-2:2] .= 10:10:50
    @time update_u!(var)
end

end

O.main()
O.main()
O.main()

# %%
@code_warntype O.main()

# %%
module P

using OffsetArrays

mutable struct Variables{U<:OffsetArray{ComplexF64, 1}}
    u::U
end

function update_u!(var::Variables)
    rhs = OffsetArray{ComplexF64}(undef, -2:2)
    rhs[-2:2] = randn(ComplexF64, 5)
    var.u .+= rhs
end

function main()
    u = OffsetArray{ComplexF64}(undef, -2:2)   
    var = Variables(u)
    var.u[-2:2] .= 10:10:50
    @time update_u!(var)
end

end

P.main()
P.main()
P.main()

# %%
@code_warntype P.main()

# %%
module Q

using OffsetArrays
using Random

struct Variables{U}
    u::U
end

function update_u!(var::Variables, rhs)
    rand!(rhs)
    var.u .+= rhs
    var
end

function main()
    u = OffsetVector{ComplexF64}(undef, -2:2)
    var = Variables(u)
    parent(var.u) .= 10:10:50
    rhs = similar(u)
    @time update_u!(var, rhs)
end

end

Q.main()
Q.main()
var = Q.main()
@show var
var.u

# %%
@code_warntype Q.main()

# %%
