height, width = 187, 746;
org_sized = rand(Float32, (2001, 2001)) * 60;
shadow_time_hrs = zeros(Float32, size(org_sized));
height_mat = rand(Float32, (height, width)) * 100; # originally values getting larger from (0, width//2) to the outside with the distance squared


angle_mat = round.(Int32, 2 .* atand.(0:height-1, (0:width-1)' .-(width/2-1))) .+ 1;

enlarged_size = zeros(eltype(org_sized), size(org_sized, 1) + height, size(org_sized, 2) + width);
enlarged_size[1:size(org_sized, 1), range(Int(width/2), length=size(org_sized, 2))] = org_sized;

weights = rand(Float32, 361)/ 10 .+ Float32(0.01);  # weights originally larger in the middle

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

function computeAllLines(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights)
    shadowed_segments = zeros(eltype(weights), 361)
    for x in 1:size(org_sized, 2) - 1
        for y  in 1:size(org_sized, 1) - 1
            shadow_time_hrs[x, y] = computeSumHours!(org_sized, enlarged_size, angle_mat, height_mat, weights, y, x, shadowed_segments)
        end
    end
    return shadow_time_hrs
end

@time result = computeAllLines(org_sized, enlarged_size, angle_mat, height_mat, shadow_time_hrs, weights);
unique(result)
