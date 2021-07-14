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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %%
InteractiveUtils.code_llvm(@nospecialize(f), @nospecialize(types=Tuple); raw=false, dump_module=false, optimize=true, debuginfo::Symbol=:default) = sprint() do io
    code_llvm(io, f, types; raw=raw, dump_module=dump_module, optimize=optimize, debuginfo=debuginfo)
end

# %%
s = @code_llvm debuginfo=:none Base.Math.sin_kernel(0.5)

# %%
s |> typeof

# %%
print(s)

# %%
