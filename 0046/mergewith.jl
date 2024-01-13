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

# %%
using BenchmarkTools

# %%
# 仮定：P,Q,R,Sのすべてのkeysは同じ。
P = Dict("a"=>1, "b"=>2, "c"=>3)
Q = Dict("a"=>2, "b"=>4, "c"=>2)
R = Dict("a"=>2, "b"=>4, "c"=>5)
S = Dict("a"=>4, "b"=>9, "c"=>1)
T = [Q, R, S]

@btime P = mergewith(+, $T...)

# %%
# 仮定：P,Q,R,Sのすべてのkeysは同じ。
P = Dict("a"=>1, "b"=>2, "c"=>3)
Q = Dict("a"=>2, "b"=>4, "c"=>2)
R = Dict("a"=>2, "b"=>4, "c"=>5)
S = Dict("a"=>4, "b"=>9, "c"=>1)
T = [Q, R, S]

for k in keys(P)
    P[k] = sum(D -> D[k], T)
end
P

# %%
# 仮定：P,Q,R,Sのすべてのkeysは同じ。
P = Dict("a"=>1, "b"=>2, "c"=>3)
Q = Dict("a"=>2, "b"=>4, "c"=>2)
R = Dict("a"=>2, "b"=>4, "c"=>5)
S = Dict("a"=>4, "b"=>9, "c"=>1)
T = [Q, R, S]

@btime begin
    for k in keys(PP)
        PP[k] = sum(D -> D[k], TT)
    end
    PP
end setup = begin
    PP = copy(P)
    TT = deepcopy(T)
end

# %%
# 仮定：PはQ,R,Sのすべてのkeysを持つ。
P = Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4, "e"=>5)
Q = Dict("a"=>2, "b"=>4, "c"=>2)
R = Dict("a"=>2, "b"=>4, "d"=>5)
S = Dict("a"=>4, "b"=>9, "e"=>1)
T = [Q, R, S]

@btime P = mergewith(+, $T...)

# %%
# 仮定：PはQ,R,Sのすべてのkeysを持つ。
P = Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4, "e"=>5)
Q = Dict("a"=>2, "b"=>4, "c"=>2)
R = Dict("a"=>2, "b"=>4, "d"=>5)
S = Dict("a"=>4, "b"=>9, "e"=>1)
T = [Q, R, S]

for k in keys(P)
    P[k] = sum(D -> get(D, k, 0), T)
end
P

# %%
#仮定：PはQ,R,Sのすべてのkeysを持つ。
P = Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4, "e"=>5)
Q = Dict("a"=>2, "b"=>4, "c"=>2)
R = Dict("a"=>2, "b"=>4, "d"=>5)
S = Dict("a"=>4, "b"=>9, "e"=>1)
T = [Q, R, S]

@btime begin
    for k in keys(PP)
        PP[k] = sum(D -> get(D, k, 0), TT)
    end
    PP
end setup = begin
    PP = copy(P)
    TT = deepcopy(T)
end

# %%
