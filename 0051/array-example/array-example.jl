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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %%
using Distributions

function alpha_level(n, C; p=1/2)
    z² = -2log(C) + log((n+1)^2/n) - log(2π*p*(1-p))
    2ccdf(Normal(), √z²)
end
strfrac(x) = isinteger(x) ? string(x) : replace(string(x), r"//"=>"/")

ps = (0.5, 0.1)
Cs = (1, 1//3, 1//10, 1//30)
ns = (10, 30, 100, 300, 1000, 3000)

for p in ps
    println("-------- p = $p")
    println("\\begin{array}{cc|$("c" ^ length(ns))}")
    println("\\multicolumn{2}{c|}{\\multirow{2}{*}{\$p_0=$p\$}}")
    println(" & \\multicolumn{$(length(ns))}{c}{n} \\\\")
    println("\\cline{3-$(2+length(ns))}")
    print(" &")
    for n in ns
        print(" & ", n)
    end
    println(" \\\\")
    println("\\hline")
    for (i, C) in enumerate(Cs)
        if i == 1
            println("\\multirow{$(length(Cs))}{*}{\$C\$}")
        end
        print(" & \\multicolumn{1}{|c|}{$(strfrac(C))}")
        for (j, n) in enumerate(ns)
            print(" & ")
            alpha_n = round(100alpha_level(n, C; p); digits=2)
            print(alpha_n, "\\%")
        end
        println(" \\\\")
    end
    println("\\end{array}")

    [round(100alpha_level(n, C; p); digits=2) for C in Cs, n in ns] |> display
    println()
end

# %%
