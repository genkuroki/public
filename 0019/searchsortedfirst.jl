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
using BenchmarkTools
using Base.Sort

# %%
function f!(a, b)
    for x in b
        push!(a, x)
        sort!(a)
    end
    a
end

function g!(a, b)
    append!(a, b)
    sort!(a)
end

function h!(a, b)
    for x in b
        i = Base.Sort.searchsortedfirst(a, x)
        insert!(a, i, x)
    end
    a
end

# %%
a, b = sort(rand(1000)), rand(10)

# %%
f!(copy(a), b) == g!(copy(a), b) == h!(copy(a), b)

# %%
@btime f!(A, $b) setup=(A = copy(a));

# %%
@btime g!(A, $b) setup=(A = copy(a));

# %%
@btime h!(A, $b) setup=(A = copy(a));

# %%
@btime f!(A, $b) setup=(A = copy(a); sizehint!(A, length(a)+length(b)));

# %%
@btime g!(A, $b) setup=(A = copy(a); sizehint!(A, length(a)+length(b)));

# %%
@btime h!(A, $b) setup=(A = copy(a); sizehint!(A, length(a)+length(b)));

# %%
A = copy(a)
@btime begin
    for i in 1:10
        push!(A, b[i])
        sort!(A)
    end
    A
end

# %%
@btime begin
    for i in 1:10
        push!(B, b[i])
        sort!(B)
    end
    B
end setup=(B = copy(a); sizehint!(B, length(a)+length(b)));

# %%
