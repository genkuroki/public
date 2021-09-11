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

# %% [markdown]
# https://discourse.julialang.org/t/how-to-read-only-the-last-line-of-a-file-txt/68005/6

# %%
str = """
some header option
header other option
random metadata
A,B,C,X,Y,Z
1,1,1,2.0,0.0,102.0
1,1,2,2.0,0.0,202.0
1,1,3,2.0,0.0,302.0
1,2,1,3.0,1.0,103.0
1,2,2,3.0,1.0,203.0
1,2,3,3.0,1.0,303.0
1,3,1,4.0,2.0,104.0
1,3,2,4.0,2.0,204.0
1,3,3,4.0,2.0,304.0
1,4,1,5.0,3.0,105.0
1,4,2,5.0,3.0,205.0
1,4,3,5.0,3.0,305.0
1,5,1,6.0,4.0,106.0
1,5,2,6.0,4.0,206.0
1,5,3,6.0,4.0,306.0
1,6,1,7.0,5.0,107.0
1,6,2,7.0,5.0,207.0
1,6,3,7.0,5.0,307.0
"""

write("test.txt", str)

# %%
function read_last(file)
    open(file) do io
        seekend(io)
        seek(io, position(io) - 1)
        while Char(peek(io)) != '\n'
            seek(io, position(io) - 1)
        end
        read(io, Char)
        read(io, String)
    end
end

# %%
function read_last2(file)
    open(file) do io
        seekend(io)
        seek(io, position(io) - 1)
        p = position(io)
        while Char(peek(io)) != '\n'
            seek(io, position(io) - 1)
        end
        if position(io) == p
            seek(io, position(io) - 1)
            while Char(peek(io)) != '\n'
                seek(io, position(io) - 1)
            end
        end
        read(io, Char)
        read(io, String)
    end |> chomp
end

# %%
read_last("test.txt")

# %%
read_last2("test.txt")

# %%
