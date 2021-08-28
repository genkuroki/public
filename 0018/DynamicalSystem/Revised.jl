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
using DynamicalSystems, DelimitedFiles, .Threads,  BenchmarkTools, ProgressMeter
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
function dds_jac(x, par, n)
    r, k = par;
    a = 5.0; mu = 0.5; d = 0.2;
    J11 = exp(- x[2]/(x[1]^2 + a) - r*(x[1]/k - 1)) - x[1]*exp(- x[2]/(x[1]^2 + a) - r*(x[1]/k - 1))*(r/k - (2*x[1]*x[2])/(x[1]^2 + a)^2)
    J12 = -(x[1]*exp(- x[2]/(x[1]^2 + a) - r*(x[1]/k - 1)))/(x[1]^2 + a)
    J21 = x[2]*exp((mu*x[1])/(x[1]^2 + a) - d)*(mu/(x[1]^2 + a) - (2*mu*x[1]^2)/(x[1]^2 + a)^2)
    J22 = exp((mu*x[1])/(x[1]^2 + a) - d)
    return @SMatrix [J11  J12; J21 J22]
end

# %%
function traj!(tr, f, p, x0, T; Ttr=0)
    x = x0
    for i in 1:Ttr
        x = f(x, p, i-1)
    end
    @inbounds tr[1] = x
    for i in 1:T
        @inbounds tr[i+1] = f(tr[i], p, Ttr+i-1)
    end
    tr
end

# %%
function traj2d!(tr1, tr2, f, p, x0, T; Ttr=0)
    x = x0
    for i in 1:Ttr
        x = f(x, p, i-1)
    end
    @inbounds tr1[1], tr2[1] = x
    for i in 1:T
        @inbounds tr1[i+1], tr2[i+1] = f((tr1[i], tr2[i]), p, Ttr+i-1)
    end
    tr1, tr2
end

# %%
traj2d!(zeros(21), zeros(21), ((a, b), p, t) -> (b, a+b), nothing, (0, 1), 20)

# %%
function seqper_new(x; tol=1e-3)       # function to calculate periodicity
    n = length(x)
    @inbounds for k in 2:(n ÷ 2 + 1)
        if abs(x[k] - x[1]) ≤ tol
            all(j -> abs(x[j] - x[j-k+1]) ≤ tol, k:n) && return k - 1
        end
    end 
    return n
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

x = 1:3
y = 1:4
@show meshh(x, y)
meshh(x, y) .== vec.(reim(complex.(x', y)))

# %%
function meshh!(xmesh, ymesh, x, y)    # create a meshgrid of x and y
    len_x, len_y = length(x), length(y)
    for i in 1:len_x
        for j in 1:len_y
            xmesh[len_y*(i-1) + j] = x[i]
            ymesh[len_y*(i-1) + j] = y[j]
        end
    end
end

x = 1:3
y = 1:4
xmesh = Vector{Float64}(undef, length(x)*length(y))
ymesh = Vector{Float64}(undef, length(x)*length(y))
meshh!(xmesh, ymesh, x, y)
@show xmesh, ymesh
meshh(x, y) .== (xmesh, ymesh)

# %%
function isoperiodic_test_org(ds, par_area::Matrix{Float64}, NIter::Int64,
                         nxblock::Int64, nyblock::Int64, NTr::Int64, xpts::Int64, ypts::Int64)    

    p1_st =  par_area[1];          # beginning of first parameter 
    p1_nd =  par_area[2];          # end of first parameter
    p2_st =  par_area[3];          # beginning of second parameter
    p2_nd =  par_area[4];          # end of second parameter

    p1_block = range(p1_st, p1_nd, length=nxblock + 1);    # divide par_area in blocks 
    p2_block = range(p2_st, p2_nd, length=nyblock + 1);    # divide par_area in blocks 
    l_bipar = xpts * ypts;                                 # total number of pts in each block
    total_blocks = nxblock * nyblock                       # total number of blocks
    sol_last = Array{Float64,2}(undef, NIter + 1, l_bipar);   # stores x values of solution for l_bipar (r,k) paris

    dss = [deepcopy(ds) for _ in 1:nthreads()]
    periods = Vector{Int}(undef, l_bipar * total_blocks)
    par1mesh = Vector{Float64}(undef, l_bipar * total_blocks)
    par2mesh = Vector{Float64}(undef, l_bipar * total_blocks)
    prog = Progress(total_blocks)
    for ii in 1:nxblock
        par1range = range(p1_block[ii], p1_block[ii+1], length=xpts+1)[1:end-1]
        for jj in 1:nyblock
            par2range = range(p2_block[jj], p2_block[jj+1], length=ypts+1)[1:end-1]
            Threads.@threads for i = 1:xpts         
                tid = threadid()
                for j = 1:ypts
                     set_parameter!(dss[tid], [par1range[i], par2range[j]]) # change the parameter values
                     tr = trajectory(dss[tid], NIter; Ttr=NTr)              # find the solution
                     sol_last[:, ypts*(i-1)+j] = tr[:,1]                    # store x-values
                end
            end
            a = l_bipar*(nyblock*(ii-1)+jj-1)
            periods[a+1:a+l_bipar] = seqper_new.(eachcol(sol_last), tol=0.001) # seqper_new calculates periodicity
            par1mesh[a+1:a+l_bipar], par2mesh[a+1:a+l_bipar] = meshh(par1range, par2range) # meshgrid of parameter values
            next!(prog)
        end
    end
    par1mesh, par2mesh, periods
end

# %%
function isoperiodic_test_rev1(ds, par_area, NIter, nxblock, nyblock, NTr, xpts, ypts)    

    p1_st = par_area[1]          # beginning of first parameter 
    p1_nd = par_area[2]          # end of first parameter
    p2_st = par_area[3]          # beginning of second parameter
    p2_nd = par_area[4]          # end of second parameter

    p1_block = range(p1_st, p1_nd, length = nxblock + 1)   # divide par_area in blocks 
    p2_block = range(p2_st, p2_nd, length = nyblock + 1)   # divide par_area in blocks 
    l_bipar = xpts * ypts                                  # total number of pts in each block
    total_blocks = nxblock * nyblock                       # total number of blocks
    sol_last = Matrix{Float64}(undef, NIter + 1, l_bipar)  # stores x values of solution for l_bipar (r,k) paris

    param = [zeros(2) for _ in 1:nthreads()]
    dss = [deepcopy(ds) for _ in 1:nthreads()]
    periods = Vector{Int}(undef, l_bipar * total_blocks)
    par1mesh = Vector{Float64}(undef, l_bipar * total_blocks)
    par2mesh = Vector{Float64}(undef, l_bipar * total_blocks)
    prog = Progress(total_blocks)
    @inbounds for ii in 1:nxblock
        par1range = range(p1_block[ii], p1_block[ii+1], length=xpts+1)[1:end-1]
        for jj in 1:nyblock
            par2range = range(p2_block[jj], p2_block[jj+1], length=ypts+1)[1:end-1]
            Threads.@threads for i = 1:xpts
                tid = threadid()
                @inbounds for j = 1:ypts
                    param[tid] .= (par1range[i], par2range[j])
                    set_parameter!(dss[tid], param[tid])                    # change the parameter values
                    tr = trajectory(dss[tid], NIter; Ttr=NTr)                # find the solution
                    sol_last[:, ypts*(i-1)+j] .= first.(tr.data) # store x-values
                end
            end
            a = l_bipar*(nyblock*(ii-1)+jj-1)
            periods[a+1:a+l_bipar] .= seqper_new.(eachcol(sol_last), tol=0.001) # seqper_new calculates periodicity
            @views meshh!(par1mesh[a+1:a+l_bipar], par2mesh[a+1:a+l_bipar], par1range, par2range) # meshgrid of parameter values
            next!(prog)
        end
    end
    par1mesh, par2mesh, periods
end

# %%
function isoperiodic_test_rev2_old(init, par_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
    u0 = SVector{2}(init)
    p1_st = par_area[1]          # beginning of first parameter 
    p1_nd = par_area[2]          # end of first parameter
    p2_st = par_area[3]          # beginning of second parameter
    p2_nd = par_area[4]          # end of second parameter

    p1_block = range(p1_st, p1_nd, length = nxblock + 1)   # divide par_area in blocks 
    p2_block = range(p2_st, p2_nd, length = nyblock + 1)   # divide par_area in blocks 
    l_bipar = xpts * ypts                                  # total number of pts in each block
    total_blocks = nxblock * nyblock                       # total number of blocks
    sol_last = Matrix{Float64}(undef, NIter + 1, l_bipar)  # stores x values of solution for l_bipar (r,k) paris

    periods = Vector{Int}(undef, l_bipar * total_blocks)
    par1mesh = Vector{Float64}(undef, l_bipar * total_blocks)
    par2mesh = Vector{Float64}(undef, l_bipar * total_blocks)
    tr = [Vector{SVector{2, Float64}}(undef, NIter + 1) for _ in 1:nthreads()]
    prog = Progress(total_blocks)
    @inbounds for ii in 1:nxblock
        par1range = range(p1_block[ii], p1_block[ii+1], length=xpts+1)[1:end-1]
        for jj in 1:nyblock
            par2range = range(p2_block[jj], p2_block[jj+1], length=ypts+1)[1:end-1]   
            Threads.@threads for i = 1:xpts
                tid = threadid()
                @inbounds for j = 1:ypts
                    traj!(tr[tid], dds_rule, (par1range[i], par2range[j]), u0, NIter; Ttr=NTr) # find the solution
                    sol_last[:, ypts*(i-1)+j] .= first.(tr[tid])   # store x-values
                end
            end
            a = l_bipar*(nyblock*(ii-1)+jj-1)
            periods[a+1:a+l_bipar] .= seqper_new.(eachcol(sol_last), tol=0.001) # seqper_new calculates periodicity
            @views meshh!(par1mesh[a+1:a+l_bipar], par2mesh[a+1:a+l_bipar], par1range, par2range) # meshgrid of parameter values
            next!(prog)
        end
    end
    par1mesh, par2mesh, periods
end

# %%
function isoperiodic_test_rev2(init, par_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
    u0 = SVector{2}(init)
    p1_st = par_area[1]          # beginning of first parameter 
    p1_nd = par_area[2]          # end of first parameter
    p2_st = par_area[3]          # beginning of second parameter
    p2_nd = par_area[4]          # end of second parameter

    p1_block = range(p1_st, p1_nd, length = nxblock + 1)   # divide par_area in blocks 
    p2_block = range(p2_st, p2_nd, length = nyblock + 1)   # divide par_area in blocks 
    l_bipar = xpts * ypts                                  # total number of pts in each block
    total_blocks = nxblock * nyblock                       # total number of blocks

    periods = Vector{Int}(undef, l_bipar * total_blocks)
    par1mesh = Vector{Float64}(undef, l_bipar * total_blocks)
    par2mesh = Vector{Float64}(undef, l_bipar * total_blocks)
    tr1 = [Vector{Float64}(undef, NIter + 1) for _ in 1:nthreads()]
    tr2 = [Vector{Float64}(undef, NIter + 1) for _ in 1:nthreads()]
    prog = Progress(total_blocks)
    @inbounds for ii in 1:nxblock
        par1range = range(p1_block[ii], p1_block[ii+1], length=xpts+1)[1:end-1]
        for jj in 1:nyblock
            par2range = range(p2_block[jj], p2_block[jj+1], length=ypts+1)[1:end-1]
            a = l_bipar*(nyblock*(ii-1)+jj-1)
            Threads.@threads for i = 1:xpts
                tid = threadid()
                @inbounds for j = 1:ypts
                    traj2d!(tr1[tid], tr2[tid], 
                        dds_rule, (par1range[i], par2range[j]), u0, NIter; Ttr=NTr) # find the solution
                    periods[a + ypts*(i-1)+j] = seqper_new(tr1[tid], tol=0.001)
                end
            end
            @views meshh!(par1mesh[a+1:a+l_bipar], par2mesh[a+1:a+l_bipar], par1range, par2range) # meshgrid of parameter values
            next!(prog)
        end 
    end
    par1mesh, par2mesh, periods
end

# %%
##  Test

parameter_area = [1.0 5.0 2.0 5.0]
nxblock = 1
nyblock = 1
xpts = 100
ypts = 100
NIter = 2000
NTr = 50000
init = [0.4, 0.5]
ds = dds_constructor(init)
##
println("********** Minor correction of the original:")
result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Revised 1:")
result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Revised 2:")
result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Confirmation of equivalence:")
@show result_org .== result_rev1;
@show result_rev1 .== result_rev2;

# %% tags=[]
##  Test

parameter_area = [1.0 5.0 2.0 5.0]
nxblock = 2
nyblock = 2
xpts = 100
ypts = 100
NIter = 2000
NTr = 50000
init = [0.4, 0.5]
ds = dds_constructor(init)
##
println("********** Minor correction of the original:")
result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Revised 1:")
result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Revised 2:")
result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Confirmation of equivalence:")
@show result_org .== result_rev1;
@show result_rev1 .== result_rev2;

# %%
##  Test

parameter_area = [1.0 5.0 2.0 5.0]
nxblock = 10
nyblock = 10
xpts = 10
ypts = 10
NIter = 2000
NTr = 50000
init = [0.4, 0.5]
ds = dds_constructor(init)
##
println("********** Minor correction of the original:")
result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
#result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
#result_org = @time isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Revised 1:")
result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
#result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
#result_rev1 = @time isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Revised 2:")
result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
#result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
#result_rev2 = @time isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)
println("********** Confirmation of equivalence:")
@show result_org .== result_rev1;
@show result_rev1 .== result_rev2;

# %%
using Plots
xoxp1(x) = x/(x+1)

r, k, period = @time isoperiodic_test_rev2([0.4, 0.5], [1.0 5.0 2.0 5.0], 2000, 1, 1, 50000, 80, 60)
r = reshape(r, 60, 80)
k = reshape(k, 60, 80)
period = reshape(period, 60, 80)
heatmap(vec(r[1,:]), vec(k[:,1]), xoxp1.(period); xlabel="r", ylabel="k", title="period/(period + 1)")

# %%
@btime isoperiodic_test_rev2($([0.4, 0.5]), $([1.0 5.0 2.0 5.0]), 2000, 1, 1, 50000, 80, 60)

# %%
@code_warntype isoperiodic_test_org(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)

# %%
@code_warntype isoperiodic_test_rev1(ds, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)

# %%
@code_warntype isoperiodic_test_rev2(init, parameter_area, NIter, nxblock, nyblock, NTr, xpts, ypts)

# %%
