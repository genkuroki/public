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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
]activate .

# %%
#using LinearAlgebra, SparseArrays, Plots, DifferentialEquations
#using Surrogates
using Plots
#gr()
L=20
N=100
k=50
x=range(0, stop=100, length=100);
dx=x[2]-x[1];
dt=0.5*dx^2/(2*k);
t=0:dt:10
# print("the vlues of t",t)
#T=zeros(length(t),length(x))
T=zeros(length(x),length(t))
s=size(T)
print("size of matrix",s)
#T[1,2]=T[1,1]=100 #complete 1st row
#T[end,end]=150  #last row of a matrix
@. T[:, 1] = exp(-(x-50)^2/(2*5^2))
for j in 1:1960
    for i=2:99
        T[i,j+1]=T[i,j]+k*dt/(dx^2)*(T[i+1,j]-2*T[i,j] + T[i-1,j])
    end
end
surface(t, x, T; color=:gist_earth)

# %%
#using LinearAlgebra, SparseArrays, Plots, DifferentialEquations
#using Surrogates
using Plots
#gr()
L=20
N=100
k=50
x=range(0, stop=100, length=100);
dx=x[2]-x[1];
dt=0.5*dx^2/(2*k);
t=0:dt:10
# print("the vlues of t",t)
#T=zeros(length(t),length(x))
T=zeros(length(x),length(t)) # t <--> x
s=size(T)
print("size of matrix",s)
#T[1,2]=T[1,1]=100 #complete 1st row
#T[end,end]=150  #last row of a matrix
@. T[:, 1] = exp(-(x-50)^2/(2*5^2)) # initial values
#for i in 2:1960
for j in 1:1960 # i -> j
#    for j=1:99
    for i=2:99 # j -> i, 1:99 -> 2:99
        T[i,j+1]=T[i,j]+k*dt/(dx^2)*(T[i+1,j]-2*T[i,j] + T[i-1,j])
    end
end
#surface(x,t,T)
surface(t,x,T) # x <--> t

# %%
T

# %%
using Plots

k = 50
x = range(0, 100; length=100);
dx = step(x)
dt = 0.5*dx^2/(2k)
t = 0:dt:10
T = zeros(length(x), length(t))
U = zeros(length(x), length(t))

@. T[:, 1] = exp(-(x-50)^2/(2*5^2))

for i in eachindex(x)
    U[i, 1] = exp(-(x[i]-50)^2/(2*5^2))
end

@show T ≈ U

plot(x, T[:, 1]; label="T[:,1]", lw=2)
plot!(x, U[:, 1]; label="U[:,1]", lw=2, ls=:dash)

# %%
@show T[:, 1];

# %%
plot(T[:, 1])

# %%
?@.

# %%
