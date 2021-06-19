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
a = quote
    x = "hoge!"
    for i in 1:5
        println(x)
    end
    println("after for loop: ", x)
end |> Base.remove_linenums! |> string
include_string(Main, a)

# %%
b = quote
    x = "hoge!"
    for i in 1:5
        x = "moge!"
        println(x)
    end
    println("after for loop: ", x)
end |> Base.remove_linenums! |> string
include_string(Main, b)

# %%
c = quote
    x = "hoge!"
    for i in 1:5
        local x = "moge!"
        println(x)
    end
    println("after for loop: ", x)
end |> Base.remove_linenums! |> string
include_string(Main, c)

# %%
d = quote
    x = "hoge!"
    for i in 1:5
        global x = "moge!"
        println(x)
    end
    println("after for loop: ", x)
end |> Base.remove_linenums! |> string
include_string(Main, d)

# %%
e = quote
    x = "hoge!"
    for i in 1:5
        x *= "moge!"
        println(x)
    end
    println("after for loop: ", x)
end |> Base.remove_linenums! |> string
include_string(Main, e)

# %%
write("b.jl", b)
read("b.jl", String) |> println

# %%
; d:/Julia-1.6.1/bin/julia b.jl

# %%
x = "hoge!"
for i in 1:5
    println(x)
end
println("after for loop: ", x)

# %%
x = "hoge!"
for i in 1:5
    x = "moge!"
    println(x)
end
println("after for loop: ", x)

# %%
for i in 1:5
    y = "moge!"
    println(x)
end
println("after for loop: ", y)

# %%
x = "hoge!"
for i in 1:5
    x *= "moge!"
    println(x)
end
println("after for loop: ", x)

# %%
