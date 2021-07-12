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
# https://discourse.julialang.org/t/using-a-abstractunitrange-inside-a-struct/64463/5

# %%
using BenchmarkTools

struct Instance_RH{P <: AbstractFloat, Q <: AbstractUnitRange}
    cost::Array{P, 2}
    time_by_quantity::Array{P, 2}
    penalty_by_multitask::Array{P, 1}
    penalty_by_know_how::Array{P, 2}
    maximum_budget::P
    I::Q
    J::Q
end

function generate_instance(cost, time_by_quantity, 
    penalty_by_multitask, penalty_by_know_how, maximum_budget)
    
    dimensions_of_instance = size(cost);
    I = Base.OneTo(dimensions_of_instance[1]);
    J = Base.OneTo(dimensions_of_instance[2]);
    return Instance_RH(cost, time_by_quantity, 
        penalty_by_multitask, penalty_by_know_how, maximum_budget, I, J);
end


cost = [404.16130257432536 342.6121400760149 334.65400017129275; 609.9826980212956 451.2171757784507 483.83148674153233; 362.93218686863423 414.41083606045623 675.6532058287046; 498.73524270040207 379.8833522296185 372.12075326006203; 507.1394146512986 343.40004673720284 344.87937542652276; 614.5935044087983 343.38652475868105 666.6455577870877; 654.5084251471844 460.7325334872372 343.51651476383165; 514.5598529059743 445.47832539448785 667.5994573488664; 352.0625777996967 357.63523731788064 426.5800063667885; 633.8737992123354 575.4684437771787 526.8544116378721]
time_by_quantity = [551.149318836251 353.5576871795756 544.7850084991596; 540.8786124298001 347.11253457830435 534.2360337745861; 530.9324781581835 340.806606793647 524.2600861079036; 521.2428271420944 334.758141874596 514.891286519538; 512.0603850942899 329.686814227325 506.0319653251741; 503.2315108697694 324.6932839541254 497.68125911747074; 495.06705853446783 320.35799835061465 489.867352509886; 487.4786234404495 316.2110977723632 483.3234557345035; 481.19585646796173 312.4731479370266 477.0680377885712; 476.0273928951479 309.1389856638628 472.01244831930615]
penalty_by_multitask = [3.932768925120902, 2.1629197821510475, 3.747201375134096, 1.2027401788064473, 3.762680990567876, 3.601216190081945, 3.9412748264981996, 1.004961550315905, 3.963908498229274, 1.5845709053894907]
penalty_by_know_how = [0.811984335664391 0.2852506669586609 0.4118216702970811; 0.10635146869948212 0.3497031430472351 0.6630721652898788; 0.13880963807944877 0.5621105391563797 0.22663788295134837; 0.44358375556142055 0.22274423735250545 0.5015070457420686; 0.14826541117764483 0.7966779051319146 0.6769308830979445; 0.07009524283219473 0.08771913744412921 0.5220504363220271; 0.8324525095520683 0.8185583447816812 0.4939598817833735; 0.597792860354574 0.4440247339167606 0.8233927217289345; 0.014589514521600778 0.29122819336562134 0.6274289868926287; 0.21229629326964108 0.2910770390826134 0.8959322834368744]
maximum_budget = 14109.108387884087
instance_0 = generate_instance(cost, time_by_quantity, penalty_by_multitask, penalty_by_know_how, maximum_budget)

function calculate!(y, x, instance) # Quantidade de pessoas i envolvidas no projeto j
    @inbounds begin
        for j in instance.J
            y[j] = 0;
            for i in instance.I
                y[j] += x[i, j];
            end
        end
    end
    return minimum(y) == 0 ? count(i->(i == 0), y) : 0;
end
x = [1 1 0; 0 0 1; 1 0 1; 0 0 0; 0 0 0; 1 1 1; 0 0 0; 0 0 0; 0 0 0; 0 1 1]
y = zeros(Int, size(x, 2))
@benchmark calculate!($y, $x, $instance_0)

# %%
function calculate_2!(y, x, instance) # Quantidade de pessoas i envolvidas no projeto j
    @inbounds begin
        for j in 1:3
            y[j] = 0;
            for i in 1:10
                y[j] += x[i, j];
            end
        end
    end
    return minimum(y) == 0 ? count(i->(i == 0), y) : 0;
end

@benchmark calculate_2!($y, $x, $instance_0)

# %%
function calculate_3!(y, x, instance) # Quantidade de pessoas i envolvidas no projeto j
    m, n = size(x)
    @inbounds begin
        for j in 1:n
            y[j] = 0;
            for i in 1:m
                y[j] += x[i, j];
            end
        end
    end
    return minimum(y) == 0 ? count(i->(i == 0), y) : 0;
end

@benchmark calculate_3!($y, $x, $instance_0)

# %%
calculate!(y, x, instance_0) == calculate_2!(y, x, instance_0) == calculate_3!(y, x, instance_0)

# %%
using BenchmarkTools

struct Instance_RH{P <: AbstractFloat, Q <: AbstractUnitRange}
    cost::Array{P, 2}
    time_by_quantity::Array{P, 2}
    penalty_by_multitask::Array{P, 1}
    penalty_by_know_how::Array{P, 2}
    maximum_budget::P
    I::Q
    J::Q
end

function generate_instance(cost, time_by_quantity, 
    penalty_by_multitask, penalty_by_know_how, maximum_budget)
    
    dimensions_of_instance = size(cost);
    I = Base.OneTo(dimensions_of_instance[1]);
    J = Base.OneTo(dimensions_of_instance[2]);
    return Instance_RH(cost, time_by_quantity, 
        penalty_by_multitask, penalty_by_know_how, maximum_budget, I, J);
end

function calculate!(y, x, instance) # Quantidade de pessoas i envolvidas no projeto j
    @inbounds begin
        for j in instance.J
            y[j] = 0;
            for i in instance.I
                y[j] += x[i, j];
            end
        end
    end
    #return minimum(y) == 0 ? count(i->(i == 0), y) : 0;
    y
end

function calculate_2!(y, x, instance) # Quantidade de pessoas i envolvidas no projeto j
    @inbounds begin
        for j in 1:1000
            y[j] = 0;
            for i in 1:3000
                y[j] += x[i, j];
            end
        end
    end
    #return minimum(y) == 0 ? count(i->(i == 0), y) : 0;
    y
end

function calculate_3!(y, x, instance) # Quantidade de pessoas i envolvidas no projeto j
    m, n = size(x)
    @inbounds begin
        for j in 1:n
            y[j] = 0;
            for i in 1:m
                y[j] += x[i, j];
            end
        end
    end
    #return minimum(y) == 0 ? count(i->(i == 0), y) : 0;
    y
end

m, n = 3000, 1000
instance_1 = generate_instance(rand(m, n), rand(m, n), rand(n), rand(m, n), rand())
x_1 = rand([0, 1], m, n)

y_1 = zeros(Int, size(x_1, 2))
y_12 = zeros(Int, size(x_1, 2))
y_13 = zeros(Int, size(x_1, 2))
@show calculate!(y_1, x_1, instance_1) == calculate_2!(y_12, x_1, instance_1) == calculate_3!(y_13, x_1, instance_1)

@btime calculate!($y_1, $x_1, $instance_1)
@btime calculate_2!($y_12, $x_1, $instance_1)
@btime calculate_3!($y_13, $x_1, $instance_1);

# %%
