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
#     display_name: Julia 1.11.3
#     language: julia
#     name: julia-1.11
# ---

# %%
using LaTeXStrings

# See https://github.com/JuliaStrings/LaTeXStrings.jl/blob/master/src/LaTeXStrings.jl#L54-L92
@eval LaTeXStrings begin
export @K_str
latexstringraw(a...) = LaTeXString(prod(string(x) for x in a))
@doc raw"""
    K"..."

is equivalent to `LaTeXString(raw"...")`, except that
`%$` can be used for interpolation.

```jldoctest
julia> K"foo"
LaTeXString("foo")

julia> K"π = %$(float(π))"
LaTeXString("π = 3.141592653589793")
```
"""
macro K_str(s::String)
    i = firstindex(s)
    buf = IOBuffer(maxsize=ncodeunits(s))
    ex = Expr(:call, GlobalRef(LaTeXStrings, :latexstringraw))
    while i <= ncodeunits(s)
        c = @inbounds s[i]
        i = nextind(s, i)
        if c === '%' && i <= ncodeunits(s)
            c = @inbounds s[i]
            if c === '$'
                position(buf) > 0 && push!(ex.args, String(take!(buf)))
                atom, i = parseatom(s, nextind(s, i), filename=string(__source__.file))
                Meta.isexpr(atom, :incomplete) && error(atom.args[1])
                atom !== nothing && push!(ex.args, atom)
                continue
            else
                print(buf, '%')
            end
        else
            print(buf, c)
        end
    end
    position(buf) > 0 && push!(ex.args, String(take!(buf)))
    return esc(ex)
end
end

@doc @K_str

# %%
L"foo"

# %%
K"foo"

# %%
L"π = %$(float(π))"

# %%
K"π = %$(float(π))"

# %%
using Plots
default(fmt=:png)
pgfplotsx()

# LaTeX のプリアンブルにフォント指定
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]

# 日本語を含むプロット
plot(sin; label=L"\sin x", title=K"じがじさんjigaji-san", xguide=K"じ", yguide=K"が")
plot!(tex_output_standalone=true)
savefig("test.tex")
plot!()

# %%
; cat test.tex

# %%
