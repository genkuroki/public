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
using Statistics
X = rand(Float32, 10^7)

x = copy(X); y = quantile!(x, [0.02, 0.98])
x = copy(X); z = [quantile!(x, 0.02), quantile!(x, 0.98)]
@show y == z;

# %%
x = copy(X); @time quantile!(x, [0.02, 0.98])
x = copy(X); @time quantile!(x, [0.02, 0.98])
x = copy(X); @time quantile!(x, [0.02, 0.98]);

# %%
x = copy(X); @time [quantile!(x, 0.02), quantile!(x, 0.98)]
x = copy(X); @time [quantile!(x, 0.02), quantile!(x, 0.98)]
x = copy(X); @time [quantile!(x, 0.02), quantile!(x, 0.98)];

# %%
# `@time`'s added version of https://github.com/JuliaLang/Statistics.jl/blob/master/src/Statistics.jl#L949-L959

function myquantile!(v::AbstractVector, p::Union{AbstractArray, Tuple{Vararg{Real}}};
                   sorted::Bool=false, alpha::Real=1., beta::Real=alpha)
    if !isempty(p)
        minp, maxp = extrema(p)
        print("Statistics._quantilesort!(v, sorted, minp, maxp):             ")
        @time Statistics._quantilesort!(v, sorted, minp, maxp)
    end
    print("map(x->Statistics._quantile(v, x, alpha=alpha, beta=beta), p):")
    return @time map(x->Statistics._quantile(v, x, alpha=alpha, beta=beta), p)
end

myquantile!(v::AbstractVector, p::Real; sorted::Bool=false, alpha::Real=1., beta::Real=alpha) = begin
    print("Statistics._quantilesort!(v, sorted, p, p):        ")
    @time Statistics._quantilesort!(v, sorted, p, p)
    print("Statistics._quantile(v, p, alpha=alpha, beta=beta):")
    @time Statistics._quantile(v, p, alpha=alpha, beta=beta)
end

println("myquantile!(x, [0.02, 0.98]):")
x = copy(X); u = myquantile!(x, [0.02, 0.98])
println("\n[myquantile!(x, 0.02), myquantile!(x, 0.98)]:")
x = copy(X); v = [myquantile!(x, 0.02), myquantile!(x, 0.98)]
println(); @show u == v;

# %%
x = copy(X); @time z = [quantile!(x, 0.499), quantile!(x, 0.501)]
x = copy(X); @time z = [quantile!(x, 0.499), quantile!(x, 0.501)]
x = copy(X); @time z = [quantile!(x, 0.499), quantile!(x, 0.501)];

# %%
x = copy(X); @time y = quantile!(x, [0.499, 0.501])
x = copy(X); @time y = quantile!(x, [0.499, 0.501])
x = copy(X); @time y = quantile!(x, [0.499, 0.501]);

# %%
@which quantile!(x, [0.02, 0.98])

# %%
@which quantile!(x, 0.02)

# %%
@which sort!(x, 1, length(x), Base.Sort.PartialQuickSort(2*10^5:98*10^5), Base.Sort.Forward)

# %%
