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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# # Binary heap and Dijkstra's algorithm
#
# * 黒木玄
# * 2020-06-26～2020-06-28, 2020-08-23, 2021-10-14
#
# 2020-06-28: [CSV.jl v0.7](https://discourse.julialang.org/t/ann-csv-jl-0-7-release/42162)以上に対応した.
#
# 2020-08-23: Julia v1.6.0-DEVを使うようにし, @inbounds を付けた. ([gist](https://gist.github.com/genkuroki/4b8dd01ccc931abe16246cdd8a0d75af))
#
# 2021-10-14: Julia v1.7.0-rc1 (2021-09-12) を使うようにし, `isless`, `rowvals`, `nzrange` を使うようにした. 警告: `isless` を使うようにしたら, `Float64` 型成分の配列のヒープソートが遅くなってしまった.

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#MyHeap型の定義" data-toc-modified-id="MyHeap型の定義-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>MyHeap型の定義</a></span></li><li><span><a href="#Heap-sort" data-toc-modified-id="Heap-sort-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Heap sort</a></span></li><li><span><a href="#Dijkstra's-algorithm" data-toc-modified-id="Dijkstra's-algorithm-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Dijkstra's algorithm</a></span><ul class="toc-item"><li><span><a href="#疎行列による重み付きグラフの取り扱い" data-toc-modified-id="疎行列による重み付きグラフの取り扱い-3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>疎行列による重み付きグラフの取り扱い</a></span></li><li><span><a href="#Dijkstra's-algorithm-の実装" data-toc-modified-id="Dijkstra's-algorithm-の実装-3.2"><span class="toc-item-num">3.2&nbsp;&nbsp;</span>Dijkstra's algorithm の実装</a></span></li></ul></li><li><span><a href="#Tokyo_Edgelist-example" data-toc-modified-id="Tokyo_Edgelist-example-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Tokyo_Edgelist example</a></span><ul class="toc-item"><li><span><a href="#read-CSV" data-toc-modified-id="read-CSV-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>read CSV</a></span></li><li><span><a href="#Benchmark" data-toc-modified-id="Benchmark-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>Benchmark</a></span></li><li><span><a href="#おまけ:-距離行列の成分をAny型にするとどうなるか" data-toc-modified-id="おまけ:-距離行列の成分をAny型にするとどうなるか-4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>おまけ: 距離行列の成分をAny型にするとどうなるか</a></span></li></ul></li></ul></div>

# %% [markdown]
# ## MyHeap型の定義

# %%
"""
    MyHeap{T}(value::T)

The heap structure realized on an array-like object value.
"""
struct MyHeap{T}
    value::T
end

parent_index(k) = k ÷ 2
left_index(k)   = 2k
right_index(k)  = 2k + 1

function heapifyup!(v, k)
    @inbounds while !isone(k)
        p = parent_index(k)
        isless(v[k], v[p]) || break
        v[k], v[p] = v[p], v[k]
        k = p
    end
end

function heapifydown!(v, k)
    N = length(v)
    @inbounds while true
        l = left_index(k)
        l > N && break
        r = right_index(k)
        if r > N
            isless(v[l], v[k]) || break
            v[k], v[l] = v[l], v[k]
            break
        end
        m =  isless(v[l], v[r]) ? l : r
        isless(v[m], v[k]) || break
        v[k], v[m] = v[m], v[k]
        k = m
    end
end

# function MyHeap(v)
#     N = length(v)
#     for k in 1:N
#         heapifyup!(v, k)
#     end
#     MyHeap{typeof(v)}(v)
# end

function MyHeap(v)
    N = length(v)
    for k in N:-1:1
        heapifydown!(v, k)
    end
    MyHeap{typeof(v)}(v)
end

MyHeap{T}() where T = MyHeap(T())
MyHeap() where T = MyHeap{Array{Float64, 1}}()

using Printf

function print_heap(v)
    N = length(v)
    iszero(N) && return
    L = floor(Int, log(2, N))
    for l in 0:L
        for i in 2^l:min(2^(l+1)-1, N)
            print(" "^((6 ÷ 2)*(2^(L-l)-1)))
            @printf("%6s", v[i])
            print(" "^((6 ÷ 2)*(2^(L-l)-1)))
        end
        println()
    end
end

print_heap(H::MyHeap) = print_heap(H.value)

Base.isempty(H::MyHeap) = isempty(H.value)
Base.length(H::MyHeap) = length(H.value)
Base.getindex(H::MyHeap, k) = H.value[k]

function Base.setindex!(H::MyHeap, val, k)
    v = H.value
    @inbounds v[k] = val
    @inbounds if !isone(k) && isless(v[k], v[parent_index(k)])
        heapifyup!(v, k)
        return
    end
    heapifydown!(v, k)
    val
end

"""
    push!(H::MyHeap, val)

Add val to H.
"""
function Base.push!(H::MyHeap, val)
    v = H.value
    N = length(v)
    push!(v, val)
    heapifyup!(v, N+1)
    H
end

"""
    popindex!(H::MyHeap, k)

Get and remove the k-th value of H.
"""
function popindex!(H::MyHeap, k)
    v = H.value
    @inbounds v[k], v[end] = v[end], v[k]
    val = pop!(v)
    heapifydown!(v, k)
    val
end

"""
    popfirst!(H::MyHeap)

Get and remove the first value of H.
"""
Base.popfirst!(H::MyHeap) = popindex!(H, 1)

# %%
# test of MyHeap

@show H = MyHeap(rand(1:100, 27))
print_heap(H)

# %%
# test of setindex!

H = MyHeap(rand(0:99, 31))
print_heap(H)
println()
@show H[31]
println()

v = copy(H.value)
v[31] = -99
print_heap(v)
println()
@show v[31]
println()

H[31] = -99
print_heap(H)

println()
println()
println()
println()

H = MyHeap(rand(0:99, 31))
print_heap(H)
println()
@show H[3]
println()

v = copy(H.value)
v[3] = 999
print_heap(v)
println()
@show v[3]
println()

H[3] = 999
print_heap(H)

# %%
# test of push!

H = MyHeap(rand(0:99, 20))
print_heap(H)
println()

w = copy(H.value)
push!(w, 5)
print_heap(w)
println()

push!(H, 5)
print_heap(H)

# %%
# test of popfirst!

H = MyHeap(rand(0:99, 12))
print_heap(H)
println()

for _ in 1:length(H)
    @show popfirst!(H)
    print_heap(H)
    println()
end

# %% [markdown]
# ## Heap sort

# %% [markdown]
# 上で定義したヒープを用いたソートと `SortingAlgorithms.jl` のヒープ・ソートを比較してみる.

# %%
function my_heap_sort(v)
    H = MyHeap(v)
    w = similar(v)
    @inbounds for i in 1:length(v)
        w[i] = popfirst!(H)
    end
    w
end

using SortingAlgorithms

v = rand(1:100, 10)
@show v
@show sort(v; alg=HeapSort)
@show my_heap_sort(v)
println()

X = rand(1:10^7, 10^7)
@time Y = sort(X; alg=HeapSort)
v = copy(X)
@time w = my_heap_sort(v)
w == Y

# %% [markdown]
# ## Dijkstra's algorithm

# %% [markdown]
# ### 疎行列による重み付きグラフの取り扱い

# %%
using SparseArrays
using LinearAlgebra

"""
    MyWeightedGraph{Tv, Ti}(W::SparseMatrixCSC{S, T})

Weighted graph structure given by the weight matrix W.
Assume that `W` is symmetric and off-diagonal.
"""
struct MyWeightedGraph{Tv, Ti}
    W::SparseMatrixCSC{Tv, Ti}
end

function MyWeightedGraph(W::AbstractMatrix{Tv}) where Tv
    W = sparse(Symmetric(sparse(W)))
    MyWeightedGraph{Tv, Int}(W)
end

const MyWGraph = MyWeightedGraph

weight_matrix(g::MyWeightedGraph) = g.W
nvertices(g::MyWeightedGraph) = size(g.W, 1)
vertices(g::MyWeightedGraph) = 1:nvertices(g)

nzptr(W::SparseArrays.AbstractSparseMatrixCSC, k) = view(rowvals(W), nzrange(W, k))
neighborhood(g::MyWeightedGraph, k) = nzptr(g.W, k)

# %% [markdown]
# 以下のセルで `SimpleWeightedGraphs.jl`, `LightGraphs.jl`, `GraphPlot` を読み込むが, 重み付きグラフのプロットのみに使用され, 最短経路の計算には使用されない.

# %%
using SimpleWeightedGraphs: SimpleWeightedGraphs
using LightGraphs: LightGraphs
using GraphPlot

function weighted_graph_plot(g; kwargs...)
    simplified_g = LightGraphs.SimpleGraph(g)
    nodelabel = 1:LightGraphs.nv(simplified_g)
    edgelabel = [LightGraphs.weights(g)[e.src, e.dst] for e in LightGraphs.edges(simplified_g)]
    gplot(simplified_g; nodelabel=nodelabel, edgelabel=edgelabel, kwargs...)
end

# %%
i = [1,  1, 2, 3, 4, 5, 6, 6,  7, 8]
j = [2,  5, 3, 4, 5, 6, 7, 9,  8, 9]
d = [1, 10, 2, 3, 4, 5, 6, 20, 7, 8]
W = sparse(i, j, d, 10, 10)
W = Symmetric(W)
W = sparse(W)
g = MyWGraph(W)
swg = SimpleWeightedGraphs.WGraph(W)

SimpleWeightedGraphs.weights(swg) |> Matrix |> display
weighted_graph_plot(swg; layout=circular_layout)

# %%
W = weight_matrix(g)
for i in 1:nvertices(g)
    for j in neighborhood(g, i)
        @show i, j, W[i, j]
    end
end

# %% [markdown]
# ### Dijkstra's algorithm の実装

# %% [markdown]
# 位置を表す key と出発点からの距離を表す val の組のヒープを扱うことによって, ヒープ中のデータの中で出発点からの距離が最小の位置が自動的に得られる.

# %%
struct KeyVal{Tv, Ti}
    key::Ti
    val::Tv
end

Base.isless(x::KeyVal, y::KeyVal) = isless(x.val, y.val)

# %%
# test of MyHeap of KeyVal

v = [KeyVal(k, rand(0:0.1:1)) for k in 1:10]
vv = v .|> x -> (x.key, x.val)
@show vv

H = MyHeap(v)
v = H.value

keys = [v[k].key for k in 1:10]
print_heap(keys)
println()
vals = [v[k].val for k in 1:10]
print_heap(vals)

H

# %% [markdown]
# Dijkstra's algorithm の解説については以下を参照せよ.
#
# * https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
# * https://www.youtube.com/results?search_query=dijkstra%27s+algorithm

# %%
"""
    DijkstraResult{S,T}(parent::S, distance::T)

The result of Dijkstra's algorithm searching the shortest path.
"""
struct DijkstraResult{S,T}
    parent::S
    distance::T
end

"""
    dijkstra(g::MyWeightedGraph{Tv, Ti}, start)

returns the result of Dijkstra's algorithm searching the shortest path from `start`.
"""
function dijkstra(g::MyWeightedGraph{Tv, Ti}, start) where {Tv, Ti}
    W = weight_matrix(g)
    N = nvertices(g)
    
    distance = fill(typemax(Tv), N)
    parent = zeros(Ti, N)
    visited = zeros(Bool, N)
    
    distance[start] = zero(Tv)
    visited[start] = true
    H = MyHeap([KeyVal(start, distance[start])])

    @inbounds while !isempty(H)
        i = popfirst!(H).key
        d_parent = distance[i]
        for j in neighborhood(g, i)
            d = d_parent + W[i, j]
            d ≥ typemax(Tv) && continue
            if !visited[j]
                visited[j] = true
                distance[j] = d
                parent[j] = i
                push!(H, KeyVal(j, d))
            elseif d < distance[j]
                distance[j] = d
                parent[j] = i
                push!(H, KeyVal(j, d))
            end
        end
    end

    DijkstraResult(parent, distance)
end

"""
    shortest_path(r::DijkstraResult, k)

returns a shortest path from `start` to `k`.
"""
function shortest_path(r::DijkstraResult, k)
    parent = r.parent
    distance = r.distance
    p = eltype(parent)[]
    distance[k] ≥ typemax(eltype(distance)) && return p
    push!(p, k)
    k = parent[k]
    while !iszero(k)
        pushfirst!(p, k)
        k = parent[k]
    end
    p
end

"""
    shortest_distance(r::DijkstraResult, k)

returns the shortest distance from `start` to `k`.
"""
shortest_distance(r::DijkstraResult, k) = r.distance[k]

# %%
i = [1,  1, 2, 3, 4, 5, 6, 6,  7, 8]
j = [2,  5, 3, 4, 5, 6, 7, 9,  8, 9]
d = [1, 10, 2, 3, 4, 5, 6, 20, 7, 8]
W = sparse(i, j, d, 10, 10)
W = Symmetric(W)
W = sparse(W)
g = MyWGraph(W)
swg = SimpleWeightedGraphs.WGraph(W)

weighted_graph_plot(swg; layout=circular_layout)

# %%
r = dijkstra(g, 1)

# %%
[(from=1, to=k, distance=shortest_distance(r, k), path=shortest_path(r, k)) for k in 1:10]

# %% [markdown]
# ## Tokyo_Edgelist example

# %% [markdown]
# ### read CSV

# %%
@time using CSV

# %%
# download("https://ndownloader.figshare.com/files/3663336", "Tokyo.zip")
# run(`unzip -o Tokyo.zip Tokyo_Edgelist.csv`)

# %%
@time csv = CSV.File("Tokyo_Edgelist.csv") # warmup
@show csv.names
@show csv.types
csv.columns

# %% [markdown]
# ### Benchmark

# %%
@time begin
    @time csv = CSV.File("Tokyo_Edgelist.csv")
    @time start_node = csv.START_NODE
    @time end_node = csv.END_NODE
    @time edge_length = csv.LENGTH
    @time D = sparse(start_node, end_node, edge_length)
    @time weighted_graph = MyWGraph(D)
    println()
    print("Dijkstra's algorithm:")
    @time dijkstra_result = dijkstra(weighted_graph, 1)
    println()
    @time dijkstra_path = shortest_path(dijkstra_result, 2539)
    println()
    @show dijkstra_path
    println()
end # warmup

# %%
@time begin
    @time csv = CSV.File("Tokyo_Edgelist.csv")
    @time start_node = csv.START_NODE
    @time end_node = csv.END_NODE
    @time edge_length = csv.LENGTH
    @time D = sparse(start_node, end_node, edge_length)
    @time weighted_graph = MyWGraph(D)
    println()
    print("Dijkstra's algorithm:")
    @time dijkstra_result = dijkstra(weighted_graph, 1)
    println()
    @time dijkstra_path = shortest_path(dijkstra_result, 2539)
    println()
    @show dijkstra_path
    println()
end

# %%
# Comparison between the above and the reki2000_route result.

reki2000_path = vec([1 110273 110421 110498 110574 110658 110706 110746 110865 110998 111102 111186 111371 111680 111808 111820 112054 112134 112217 112452 112503 112530 112546 112774 113010 113100 113710 113830 113926 113999 114010 114092 114308 114458 114704 114804 114969 115015 115091 115303 115352 115645 115667 115733 115771 115859 116007 116227 116606 116605 116815 117103 117447 117610 117607 117860 118052 426 510 567 617 630 666 749 833 939 956 1029 1127 1214 1350 1483 1521 1678 1742 1752 1813 1960 2013 2119 2232 2351 2408 2434 2464 2469 2539])

dijkstra_path == reki2000_path

# %% [markdown]
# See https://qiita.com/reki2000/items/55ef54b96b26d80ad694 for the reki2000_route result.

# %%
[(from=1, to=k, distance=dijkstra_result.distance[k], 
        path=shortest_path(dijkstra_result, k)) for k in 1:20]

# %%
D

# %% [markdown]
# ### おまけ: 距離行列の成分をAny型にするとどうなるか

# %%
using BenchmarkTools

Base.typemax(::Type{Any}) = typemax(Float64)
Base.zero(::Type{Any}) = zero(Float64)

D = sparse(start_node, end_node, edge_length)
D_Any = SparseMatrixCSC{Any, Int}(D)
weighted_graph = MyWGraph(D)
weighted_graph_Any = MyWGraph(D_Any)
dijkstra(weighted_graph_Any, 1)

@show typeof(weighted_graph)
@show typeof(weighted_graph_Any)
println()

print("Float64:")
@btime dijkstra(weighted_graph, 1)
print("Any:    ")
@btime dijkstra(weighted_graph_Any, 1)

# %%
D_Any

# %%
