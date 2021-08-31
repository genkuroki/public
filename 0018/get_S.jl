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
W = trues(10000, 10000)
W[rand(1:length(W), round(Int, 0.1length(W)))] .= 0
W

# %%
function get_S(W)
    zeroidxs = getindex.(findall(iszero, W), [1 2])
    S1 = Set(zeroidxs[:,1])
    S2 = Set(zeroidxs[:,2])
    return S1, S2
end

# %%
@time get_S(W)
@time get_S(W)
@time get_S(W)

# %%
function get_S1(W)
    zeroidxs = getindex.(findall(iszero, W), [1 2])
    S1 = unique(zeroidxs[:,1])
    S2 = unique(zeroidxs[:,2])
    return S1, S2
end

# %%
@time get_S1(W)
@time get_S1(W)
@time get_S1(W)

# %%
function get_S2(W)
    m, n = size(W)
    b1, b2 = falses(m), falses(n)
    for j in 1:n
        for i in 1:n
            @inbounds if W[i, j] == 0
                b1[i] = b2[j] = true
            end
        end
    end
    findall(b1), findall(b2)
end

# %%
@time get_S2(W)
@time get_S2(W)
@time get_S2(W)

# %%
sort.(get_S1(W)) == get_S2(W)

# %%
