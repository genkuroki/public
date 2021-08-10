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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
module My

struct Hello{T}
    a::T
end

end

My.Hello("Julia")

# %%
module My

struct Hello{T}
    a::T
end

function Base.show(io::IO, x::Hello)
    print(io, "Hello, ")
    show(io, x.a)
    print(io, '!')
end

end

My.Hello("Julia")

# %%
dump(My.Hello("Julia"))

# %%
π

# %%
My.Hello(π)

# %%
module My

struct Hello{T}
    a::T
end

function Base.show(io::IO, x::Hello)
    print(io, "Hello, ")
    show(io, x.a)
    print(io, '!')
end

function Base.show(io::IO, ::MIME"text/plain", x::Hello)
    print(io, "Hello, ")
    show(io, MIME"text/plain"(), x.a)
    print(io, "!")
end

end

My.Hello("Julia")

# %%
My.Hello(π)

# %%
