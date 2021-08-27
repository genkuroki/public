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

# %% [markdown]
# https://discourse.julialang.org/t/solving-difference-equation-part-2/67057

# %%
using DynamicalSystems, DelimitedFiles, .Threads,  BenchmarkTools
## Components of a test DiscreteDynamicalSystem

function dds_constructor(u0 = [0.5, 0.7];  r=1.0, k=2.0)
    return DiscreteDynamicalSystem(dds_rule, u0, [r, k], dds_jac)
end
## equations of motion:
function dds_rule(x, par, n)
   r, k = par
   a, mu, d = 5.0, 0.5, 0.2
   dx = x[1]*exp(r*(1-x[1]/k)-x[2]/(a+x[1]^2))
   dy = x[2]*exp(mu*x[1]/(a+x[1]^2)-d)
   return @SVector [dx, dy]
end
## Jacobian:
function dds_jac( x, par, n)
    r, k = par;
    a = 5.0; mu = 0.5; d = 0.2;
    J11 = exp(- x[2]/(x[1]^2 + a) - r*(x[1]/k - 1)) - x[1]*exp(- x[2]/(x[1]^2 + a) - r*(x[1]/k - 1))*(r/k - (2*x[1]*x[2])/(x[1]^2 + a)^2)
    J12 = -(x[1]*exp(- x[2]/(x[1]^2 + a) - r*(x[1]/k - 1)))/(x[1]^2 + a)
    J21 = x[2]*exp((mu*x[1])/(x[1]^2 + a) - d)*(mu/(x[1]^2 + a) - (2*mu*x[1]^2)/(x[1]^2 + a)^2)
    J22 = exp((mu*x[1])/(x[1]^2 + a) - d)
    return @SMatrix [J11  J12; J21 J22]
end

# %%
function seqper_new(x; tol=1e-3)       # function to calculate periodicity 
    @inbounds for k in 2:length(x)
        if abs(x[k] - x[1]) ≤ tol
            all(j -> abs(x[j] - x[j-k+1]) ≤ tol, k:lastindex(x)) && return k - 1
        end
    end 
    return length(x)
end

# %%
function meshh(x,y)    # create a meshgrid of x and y
    len_x = length(x);
    len_y = length(y);
    xmesh = Float64[]
    ymesh = Float64[]
    for i in 1:len_x
        for j in 1:len_y
            push!(xmesh, x[i])
            push!(ymesh, y[j])
        end
    end
    return  xmesh, ymesh 
end

# %%
#function isoperiodic_test(ds, par_area::Matrix{Float64}, NIter::Int64, filename::String,
function isoperiodic_test(ds, par_area::Matrix{Float64}, NIter::Int64, iso_filename::String,
                         nxblock::Int64, nyblock::Int64, NTr::Int64, xpts::Int64, ypts::Int64)    

    io = open(iso_filename, "a")   # opeing file to write output
    p1_st =  par_area[1];          # beginning of first parameter 
    p1_nd =  par_area[2];          # end of first parameter
    p2_st =  par_area[3];          # beginning of second parameter
    p2_nd =  par_area[4];          # end of second parameter

    p1_block = range(p1_st, p1_nd, length=nxblock + 1);    # divide par_area in blocks 
    p2_block = range(p2_st, p2_nd, length=nyblock + 1);    # divide par_area in blocks 
    l_bipar = xpts * ypts;                                 # total number of pts in each block
    total_blocks = nxblock * nyblock                       # total number of blocks
    sol_last = Array{Float64,2}(undef, NIter + 1, l_bipar);   # stores x values of solution for l_bipar (r,k) paris

    for ii in 1:nxblock
         step1size = (p1_block[ii + 1] - p1_block[ii]) / xpts;
         par1range = range(p1_block[ii], p1_block[ii + 1] - step1size, length=xpts);
        for jj in 1:nyblock
            #display("Progress: " * string(nxblock * (ii - 1) + jj) * " out of " * string(total_blocks) * " steps.")
            display("Progress: " * string(nyblock * (ii - 1) + jj) * " out of " * string(total_blocks) * " steps.")
            step2size = (p2_block[jj + 1] - p2_block[jj]) / ypts;
            par2range = range(p2_block[jj], p2_block[jj + 1] - step2size, length=ypts);   
            Threads.@threads for i = 1:xpts         
                for j = 1:ypts
                     set_parameter!(ds, [par1range[i], par2range[j]]) # change the parameter values <-- thread unsafe! bug!
                     tr = trajectory(ds, NIter; Ttr=NTr);             # find the solution
                     #sol_last[:,(i - 1) * xpts + j] =  tr[:,1];       # store x-values
                     sol_last[:,(i - 1) * ypts + j] =  tr[:,1];       # store x-values
                end
            end
            periods = seqper_new.(eachcol(sol_last), tol=0.001)       # seqper_new calculates periodicity 
            par1mesh, par2mesh = meshh(par1range, par2range)          # meshgrid of parameter values
            writedlm(io, [par1mesh par2mesh periods])                 # write data in txt file
        end 
    end
    close(io)
    return iso_filename
end

# %%
##  Test

parameter_area = [1.0 5.0 2.0 5.0]
nxblock = 1;
nyblock = 1;
xpts = 100;
ypts = 100;
NIter = 2000;
NTr = 50000;
init = [0.4, 0.5];
ds = dds_constructor(init);
##
iso_filename = "isoperiodic_test_data.txt"
isfile(iso_filename) && rm(iso_filename)
@time isoperiodic_test(ds, parameter_area, NIter, iso_filename, nxblock, nyblock, NTr, xpts, ypts)

# %%
##  Test

parameter_area = [1.0 5.0 2.0 5.0]
nxblock = 2;
nyblock = 2;
xpts = 100;
ypts = 100;
NIter = 2000;
NTr = 50000;
init = [0.4, 0.5];
ds = dds_constructor(init);
##
iso_filename = "isoperiodic_test_data2x2.txt"
isfile(iso_filename) && rm(iso_filename)
@time isoperiodic_test(ds, parameter_area, NIter, iso_filename, nxblock, nyblock, NTr, xpts, ypts)

# %%
