# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %%
function triple_dice(N::Int)
    count = 0
    for _ in 1:N
        k = rand(1:6) + rand(1:6) + rand(1:6)
        if k % 2 == 1
            count += 1
        end
    end
    return count / N
end

@time triple_dice(10^8)
@time triple_dice(10^8)
@time triple_dice(10^8)

# %%
@show N = Int32(10)^8
triple_dice(N)

# %%
function triple_dice_rev1(N::Integer)
    count = 0
    for _ in 1:N
        k = rand(1:6) + rand(1:6) + rand(1:6)
        count += k % 2 == 1
    end
    count / N
end

@time triple_dice_rev1(10^8)
@time triple_dice_rev1(10^8)
@time triple_dice_rev1(10^8)

# %%
@show N = Int32(10)^8
triple_dice_rev1(N)

# %%
