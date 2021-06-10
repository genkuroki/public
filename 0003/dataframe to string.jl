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
-8.06,-17.51,-3.55,-2.31,-10.52
5.47,-3.17,-16.43,-9.26,-3.25
-0.39,-7.39,-11.55,5.53,-1.24
11.59,1.36,-9.81,0.15,23.74
-3.41,11.27,-0.65,4.88,-10.86
-2.4,9.17,-6.35,14.36,2.94
14.98,-8.69,-6.85,1.23,-14.48
""")

# %%
df2 == df

# %%
