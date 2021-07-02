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
VERSION

# %%
using RCall

x = 1

f(var::String) = eval(Meta.parse("@rput ("*var*")"))

module M1
using RCall
import RCall: @rput 

f(var::String) = eval(Meta.parse("@rput ("*var*")"))

end

# %%
f("x")

# %%
M1.f("x")

# %%
M1.f("Main.x")

# %%
using RCall

x = 1

f(var::String) = eval(:(@rput($(Symbol(var)))))

module M2

using RCall

x = 2

f(var::String) = eval(:(@rput($(Symbol(var)))))

g(var::String) = Main.eval(:(@rput($(Symbol(var)))))
g(m::Module, var::String) = Core.eval(m, :(@rput($(Symbol(var)))))

end

# %%
M2.f("x")

# %%
R"x"

# %%
M2.g("x")

# %%
R"x"

# %%
M2.g(M2, "x")

# %%
R"x"

# %% [markdown]
# ```julia
# using RCall
#
# x = 1
#
# f(var::String) = eval(:(@rput($(Symbol(var)))))
#
# module M2
#
# using RCall
#
# x = 2
#
# f(var::String) = eval(:(@rput($(Symbol(var)))))
#
# g(var::String) = Main.eval(:(@rput($(Symbol(var)))))
#
# g(m::Module, var::String) = Core.eval(m, :(@rput($(Symbol(var)))))
#
# end
# ```
#
# ```julia
# julia> M2.f("x")
# 2
#
# julia> R"x"
# RObject{IntSxp}
# [1] 2
# ```
#
# ```julia
# julia> M2.g("x")
# 1
#
# julia> R"x"
# RObject{IntSxp}
# [1] 1
# ```
#
# ```julia
# julia> M2.g(M2, "x") # equivalent to M2.f("x")
# 2
#
# julia> R"x"
# RObject{IntSxp}
# [1] 2
# ```
#

# %%
using RCall

x = 1

f(var::String) = eval(:(@rput($(Symbol(var)))))

module M2

using RCall

x = 2

f(var::String) = eval(:(@rput($(Symbol(var)))))

g(var::String) = Main.eval(:(@rput($(Symbol(var)))))

g(m::Module, var::String) = Core.eval(m, :(@rput($(Symbol(var)))))

end

# %%
using RCall

x = 1

f(var::String) = eval(:(@rput($(Symbol(var)))))

module M2

using RCall

x = 2

f(var::String) = eval(:(@rput($(Symbol(var)))))
g(var::String) = Main.eval(:(@rput($(Symbol(var)))))
g(m::Module, var::String) = Core.eval(m, :(@rput($(Symbol(var)))))

end

# %%
