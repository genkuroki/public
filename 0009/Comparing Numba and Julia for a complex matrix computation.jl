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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://discourse.julialang.org/t/comparing-numba-and-julia-for-a-complex-matrix-computation/63703/11

# %%
versioninfo()

# %%
# resized original version

height, width = 187, 746;
#org_sized = rand(Float32, (2001, 2001)) * 60;
org_sized = rand(Float32, (401, 401)) * 60;
#org_sized = rand(Float32, (101, 101)) * 60;
shadow_time_hrs = zeros(Float32, size(org_sized));
height_mat = rand(Float32, (height, width)) * 100; # originally values getting larger from (0, width//2) to the outside with the distance squared


angle_mat = round.(Int32, 2 .* atand.(0:height-1, (0:width-1)' .-(width/2-1))) .+ 1;

enlarged_size = zeros(eltype(org_sized), size(org_sized, 1) + height, size(org_sized, 2) + width);
enlarged_size[1:size(org_sized, 1), range(Int(width/2), length=size(org_sized, 2))] = org_sized;

weights = rand(Float32, 361)/ 10 .+ 0.01;  # weights originally larger in the middle

function computeSumHours(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x)
    height, width = size(height_mat)
    short_elevations = @view enlarged_size[y:y+height, x:x+width]
    shadowed_segments = zeros(eltype(weights), 361)

    for x2 in 1:width
        for y2 in 1:height
            overshadowed = (short_elevations[y2, x2] - org_sized[y, x]) > height_mat[y2, x2]
            if overshadowed
                angle = angle_mat[y2, x2]
                if shadowed_segments[angle] == 0.0
                    shadowed_segments[angle] = weights[angle]
                end
            end
        end
    end

    return sum(shadowed_segments)
end

function computeAllLines(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights)
    for x in 1:size(org_sized, 2) - 1
        for y  in 1:size(org_sized, 1) - 1
            shadow_time_hrs[x, y] = computeSumHours(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x)
        end
    end
    return shadow_time_hrs
end

@time result = computeAllLines(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights)
@time result = computeAllLines(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights)
@time result = computeAllLines(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights)
r0 = result;

# %%
# revised version

weights = Float32.(weights)
shadowed_segments = zeros(eltype(weights), 361)

function computeSumHours!(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x, shadowed_segments)
    height, width = size(height_mat)
    short_elevations = @view enlarged_size[y:y+height, x:x+width]
    shadowed_segments .= 0

    @inbounds for x2 in 1:width
        for y2 in 1:height
            overshadowed = (short_elevations[y2, x2] - org_sized[y, x]) > height_mat[y2, x2]
            if overshadowed
                angle = angle_mat[y2, x2]
                if shadowed_segments[angle] == 0.0
                    shadowed_segments[angle] = weights[angle]
                end
            end
        end
    end

    return sum(shadowed_segments)
end

function computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
    for x in 1:size(org_sized, 2) - 1
        for y  in 1:size(org_sized, 1) - 1
            shadow_time_hrs[x, y] = computeSumHours!(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x, shadowed_segments)
        end
    end
    return shadow_time_hrs
end

@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
r1 = result;

# %%
# revised version 2

weights = Float32.(weights)
shadowed_segments = zeros(eltype(weights), 361)

function computeSumHours!(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x, shadowed_segments)
    height, width = size(height_mat)
    short_elevations = @view enlarged_size[y:y+height, x:x+width]
    shadowed_segments .= 0

    Threads.@threads for x2 in 1:width
        for y2 in 1:height
            @inbounds overshadowed = (short_elevations[y2, x2] - org_sized[y, x]) > height_mat[y2, x2]
            @inbounds if overshadowed
                angle = angle_mat[y2, x2]
                if shadowed_segments[angle] == 0.0
                    shadowed_segments[angle] = weights[angle]
                end
            end
        end
    end

    return sum(shadowed_segments)
end

function computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
    for x in 1:size(org_sized, 2) - 1
        for y  in 1:size(org_sized, 1) - 1
            shadow_time_hrs[x, y] = computeSumHours!(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x, shadowed_segments)
        end
    end
    return shadow_time_hrs
end

@show Threads.nthreads()
@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
r2 = result;

# %%
# revised version 3

weights = Float32.(weights)
shadowed_segments = zeros(eltype(weights), 361)

function computeSumHours!(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x, shadowed_segments)
    height, width = size(height_mat)
    short_elevations = @view enlarged_size[y:y+height, x:x+width]
    shadowed_segments .= 0

    Threads.@threads for x2 in 1:width
        @inbounds @simd for y2 in 1:height
            overshadowed = (short_elevations[y2, x2] - org_sized[y, x]) > height_mat[y2, x2]
            if overshadowed
                angle = angle_mat[y2, x2]
                if shadowed_segments[angle] == 0.0
                    shadowed_segments[angle] = weights[angle]
                end
            end
        end
    end

    return sum(shadowed_segments)
end

function computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
    shadowed_segments = zeros(eltype(weights), 361)
    for x in 1:size(org_sized, 2) - 1
        for y  in 1:size(org_sized, 1) - 1
            shadow_time_hrs[x, y] = computeSumHours!(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x, shadowed_segments)
        end
    end
    return shadow_time_hrs
end

@show Threads.nthreads()
@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
@time result = computeAllLines!(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights, shadowed_segments)
r3 = result;

# %%
r0 ≈ r1 ≈ r2 ≈ r3

# %%
