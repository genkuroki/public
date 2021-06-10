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
#     display_name: Julia 1.7.0-DEV
#     language: julia
#     name: julia-1.7
# ---

# %%
using DataFrames, CSV
csvstring(df) = String(take!(CSV.write(IOBuffer(), df)))
csv2df(str) = CSV.read(IOBuffer(str), DataFrame)

# %%
using Random
Random.seed!(4649373)

# %%
data = round.(10randn(7, 5), digits=2)
names = string.('a':'e')
df = DataFrame(data, names)

# %%
str = csvstring(df)

# %%
df1 = csv2df(str)

# %%
df1 == df

# %%
print(str)

# %%
# copy and paste the above

df2 = csv2df("""
a,b,c,d,e
-5.91,19.26,19.44,1.03,1.85
-3.67,-10.33,1.57,-17.67,-0.47
6.72,8.81,-3.87,-11.0,-0.56
-6.2,10.66,11.52,6.19,-7.78
-8.64,-1.6,9.28,-1.8,3.72
-1.01,-17.56,7.02,1.08,18.56
-20.48,-9.66,5.31,-17.65,11.51
""")

# %%
df2 == df

# %%
