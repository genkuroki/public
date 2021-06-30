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
using MetaUtils

# %% [markdown]
# * https://github.com/mitmath/18S096/blob/master/lectures/other/Syntax%20trees%20in%20Julia.ipynb

# %%
# https://github.com/mitmath/18S096/blob/master/lectures/other/Syntax%20trees%20in%20Julia.ipynb

using LightGraphs
using TikzGraphs
using MacroTools

"Return the number of the new vertex"
add_numbered_vertex!(g) = (add_vertex!(g); top = nv(g))

"Convert the current node into a label"
function label(sym) 
    sym == :(^) && return "\\textasciicircum"  # TikzGraphs chokes on ^
    return string("\\texttt{", sym, "}")
end

"""
Traverse the expression and build up the graph g and the set of labels for each node. (These both get modified in place.) `show_call` specifies whether to include `call` nodes in the graph in function calls. Including them represents the Julia AST more precisely, but adds a lot of visual noise to the display.

traverse! returns the number of the vertex at the top of the subexpression.
"""
function traverse!(g, labels, ex, show_call=true)
    top_vertex = add_numbered_vertex!(g)
    start_argument = 1  # which argument to start with

    if !(show_call) && ex.head == :call 
        f = ex.args[1]   # the function name
        push!(labels, label(f))
        start_argument = 2   # drop "call" from tree
    else
        push!(labels, label(ex.head))
    end
    
    for i in start_argument:length(ex.args)
        if isa(ex.args[i], Expr)
            child = traverse!(g, labels, ex.args[i], show_call)
            add_edge!(g, top_vertex, child)
        else
            n = add_numbered_vertex!(g)
            add_edge!(g, top_vertex, n)
            push!(labels, label(ex.args[i]))
        end
    end
    
    return top_vertex
end

function make_graph(ex::Expr, show_call=false)
    # https://github.com/FluxML/MacroTools.jl/blob/master/src/utils.jl#L87
    ex = MacroTools.striplines(ex)
    g = Graph()
    labels = String[]
    traverse!(g, labels, ex, show_call)
    return g, labels
end

function draw_syntax_tree(ex::Expr, show_call=false)
    TikzGraphs.plot(make_graph(ex::Expr, show_call)...)
end

macro tree(ex::Expr)
    draw_syntax_tree(ex)
end

macro tree_with_call(ex::Expr)
    draw_syntax_tree(ex, true)
end

# %%
t = @tree z = 3x + (2y+1)^2

# %%
t.data

# %%
@tree_with_call z = 3x + (2y+1)^2

# %%
@show_tree z = 3x + (2y+1)^2

# %%
@show_texpr z = 3x + (2y+1)^2

# %%
@tree z = f(g(f(x, y)))

# %%
@tree begin
    a = b + c
    d = 2a^2 + a
end

# %%
@tree (f(x^2 + 2x), g(y - (x^2 + 2x)))

# %%
@tree @time sin(10)

# %%
@show_tree @time sin(10)

# %%
@show_texpr @time sin(10)

# %%
g, labels = make_graph(:(x^2 + y^2))

for i in 1:length(labels)
    labels[i] = string("\$v_", i, "\$: \$", labels[i], "\$")
end

TikzGraphs.plot(g, labels)

# %%
@tree x*y - sin(z)

# %%
@tree (c = 0; for k in 1:n c += rand()^2 + rand()^2 <= 1 end; 4c/n)

# %%
@show_tree (c = 0; for k in 1:n c += rand()^2 + rand()^2 ≤ 1 end; 4c/n)

# %%
@show_texpr (c = 0; for k in 1:n c += rand()^2 + rand()^2 ≤ 1 end; 4c/n)

# %%
