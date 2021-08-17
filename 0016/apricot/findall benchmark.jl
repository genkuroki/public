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
using Random

n = 2^19
bitmask = bitrand(n)
boolmask = Vector{Bool}(bitmask)
bytemask = Vector{Int8}(bitmask)

@show summary(bitmask)
@show summary(boolmask)
@show summary(bytemask)
@show findall(bitmask) == findall(.!(.!bitmask)) == findall(!, .!bitmask) ==
      findall(boolmask) == findall(Bool.(boolmask)) == findall(.!(.!boolmask)) == findall(!, .!boolmask) ==
      findall(Bool.(bytemask)) == findall(bytemask .== 1) == findall(==(1), bytemask)

println(); flush(stdout)

print("findall(\$bitmask):              ")
@btime findall($bitmask)
print("findall(.!\$(.!bitmask)):        ")
@btime findall(.!$(.!bitmask))
print("findall(!, \$(.!bitmask)):     ")
@btime findall(!, $(.!bitmask))

print("findall(\$boolmask):           ")
@btime findall($boolmask)
print("findall(Bool.(\$boolmask)):      ")
@btime findall(Bool.($boolmask))
print("findall(.!\$(.!boolmask)):       ")
@btime findall(.!$(.!boolmask))
print("findall(!, \$(.!boolmask)):    ")
@btime findall(!, $(.!boolmask))

print("findall(Bool.(\$bytemask)):    ")
@btime findall(Bool.($bytemask))
print("findall(\$bytemask .== 1):       ")
@btime findall($bytemask .== 1)
print("findall(==(1), \$bytemask):    ")
@btime findall(==(1), $bytemask)
;

# %%
