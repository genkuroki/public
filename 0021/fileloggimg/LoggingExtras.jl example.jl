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
# * https://discourse.julialang.org/t/how-to-save-logging-output-to-a-log-file/14004/11
# * https://github.com/JuliaLogging/LoggingExtras.jl

# %%
#using Logging
using LoggingExtras

# %%
io = open("a.txt", "w+")
logger = SimpleLogger(io)
with_logger(logger) do
    @info(" here is some context specific logging with SimpleLogger")
end
read("a.txt", String) |> print

# %%
flush(io)
read("a.txt", String) |> print

# %%
close(io)

# %%
io = open("b.txt", "w+")
logger = FileLogger(io)
with_logger(logger) do
    @info(" here is some context specific logging with FileLogger")
end
read("a.txt", String) |> print

# %%
close(io)

# %%
?FileLogger

# %%
