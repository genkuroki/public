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
module O

struct Program f::Function end
const TestSet = Dict{Symbol, Program}()

filenames = ("dir/A.jl", "dir/B.jl", "dir/C.jl")
rmext(x) = replace(x, r"\.[^.]*$"=>"")

for fn in filenames
    s = rmext(basename(fn))
    m = Symbol(s)
    @eval module $m
        using ..O
        f() = "f defined in the module O." * $s # or include(fn)
        O.TestSet[Symbol($s)] = O.Program(f)
    end
end

end

# %%
O.TestSet

# %%
O.TestSet[:A].f()

# %%
O.A.f()

# %%
O.TestSet[:B].f()

# %%
O.TestSet[:A].f |> dump

# %%
O.TestSet[:B].f |> dump

# %%
O.TestSet[:C].f |> dump

# %%
module P

struct Program f::Function end
const TestSet = Dict{Symbol, Program}()

filenames = ("dir/A.jl", "dir/B.jl", "dir/C.jl")
rmext(x) = replace(x, r"\.[^.]*$"=>"")

for fn in filenames
    let m = Symbol(rmext(basename(fn)))
        f() = "f defined in " * string(m) * ".jl" # or include(fn)
        TestSet[m] = Program(f)
    end
end

end

# %%
P.TestSet

# %%
P.TestSet[:A].f()

# %%
P.TestSet[:B].f()

# %%
P.TestSet[:A].f |> dump

# %%
P.TestSet[:B].f |> dump

# %%
P.TestSet[:C].f |> dump

# %%
