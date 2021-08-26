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
# https://projecteuler.net/index.php?section=problems&id=30
#
# https://discourse.julialang.org/t/reduce-number-of-allocations/66990

# %%
using BenchmarkTools

# %%
function power_digit_sum(pow, n)
    return sum(c^pow for c in reverse(digits(n)))
end

function Problem30()
    ans = sum(i for i in 2:1_000_000 if i == power_digit_sum(5, i))
    return ans
end

@show Problem30();
@btime Problem30();

# %% [markdown]
# https://docs.julialang.org/en/v1/base/numbers/#Base.digits!

# %%
module Rev1

function power_digit_sum!(pow, n, ws)
    return sum(c^pow for c in digits!(ws, n))
end

function Problem30(N, m)
    ws = Vector{Int}(undef, N)
    ans = sum(i for i in 2:(10^N - 1) if i == power_digit_sum!(m, i, ws))
    return ans
end

end

@show Rev1.Problem30(5, 4);
@show Rev1.Problem30(6, 5);
@btime Rev1.Problem30(6, 5);

# %%
module DNF1

function powsum(pow, n)
    s = 0
    while n > 0
        (n, r) = divrem(n, 10)
        s += r^pow
    end
    return s
end

problem30() = sum(i for i in 2:1_000_000 if i == powsum(5, i))

end

@show DNF1.problem30();
@btime DNF1.problem30();

# %%
module DNF2

function pow5(x)
    y = x^2
    return y^2 * x
end

function pow5sum(n)
    s = 0
    while n >= 10
        (n, r) = divrem(n, 10)
        s += pow5(r)
    end
    return s + pow5(n)
end

problem30() = sum(i for i in 2:1_000_000 if i == pow5sum(i))

end

@show DNF2.problem30();
@btime DNF2.problem30();

# %% [markdown]
# https://docs.julialang.org/en/v1/devdocs/cartesian/

# %%
using Base.Cartesian

# %%
(@macroexpand @nexprs 4 d -> s += 10^(d-1)*i_d) |> Base.remove_linenums!

# %%
(@macroexpand @nloops 4 i d -> 0:9 begin
    @nexprs 4 d -> s += 10^(d-1)*i_d
end) |> Base.remove_linenums!

# %%
N = 4
s = 0
@eval @nloops $N i d -> 0:9 begin
    @nexprs $N d -> global s += 10^(d-1)*i_d
end
s

# %%
sum(0:9999)

# %%
begin
    begin
        a = Int[]
        s = 0
        x = 0
        @nloops 5 i (d -> 0:9) begin
            y = 0
            @nexprs 5 d -> y += i_d ^ 4
            if x > 1 && x == y
                push!(a, x)
                s += x
            end
            x += 1
        end
        a, s
    end
end

# %%
(@macroexpand begin
    begin
        p = @ntuple 10 i -> (i - 1) ^ 5
        a = Int[]
        s = 0
        x = 0
        @inbounds @nloops 6 i d -> 0:9 d -> y_d = p[i_d+1] begin
            y = 0
            @nexprs 6 d -> y += y_d
            if x > 1 && x == y
                push!(a, x)
                s += x
            end
            x += 1
        end
        a, s
    end
end) |> Base.remove_linenums!

# %%
@generated function _projecteuler30v(::Val{N}, ::Val{m}) where {N, m}
    quote
        p = @ntuple 10 i -> (i - 1) ^ $m
        a = Int[]
        s = 0
        x = 0
        @inbounds @nloops $N i d -> 0:9 d -> y_d = p[i_d+1] begin
            y = 0
            @nexprs $N d -> y += y_d
            if x > 1 && x == y
                push!(a, x)
                s += x
            end
            x += 1
        end
        a, s
    end
end
projecteuler30v(N, m) = _projecteuler30v(Val(N), Val(m))

@show projecteuler30v(5, 4);
@show projecteuler30v(6, 5);
@btime projecteuler30v(6, 5);

# %%
@generated function _projecteuler30(::Val{N}, ::Val{m}) where {N, m}
    quote
        p = @ntuple 10 i -> (i - 1) ^ $m
        s = 0
        x = 0
        @inbounds @nloops $N i d -> 0:9 begin
            y = 0
            @nexprs $N d -> y += p[i_d+1]
            s += (x == y) * x
            x += 1
        end
        s - 1
    end
end
projecteuler30(N, m) = _projecteuler30(Val(N), Val(m))

@show projecteuler30(5, 4);
@show projecteuler30(6, 5);
@btime projecteuler30(6, 5);

# %%
@generated function _projecteuler30x(::Val{N}, ::Val{m}) where {N, m}
    quote
        p = @ntuple 10 i -> (i - 1) ^ $m
        s = 0
        x = 0
        @inbounds @nloops $N i d -> 0:9 d -> begin
            y_d = p[i_d+1]
        end begin
            y = 0
            @nexprs $N d -> y += y_d
            s += (x == y) * x
            x += 1
        end
        s - 1
    end
end
projecteuler30x(N, m) = _projecteuler30x(Val(N), Val(m))

@show projecteuler30x(5, 4);
@show projecteuler30x(6, 5);
@btime projecteuler30x(6, 5);

# %%
function projecteuler30p()
    p = @ntuple 10 i -> (i - 1) ^ 5
    s = 0
    x = 0
    @inbounds @nloops 6 i d -> 0:9 begin
        y = 0
        @nexprs 6 d -> y += p[i_d+1]
        s += (x == y) * x
        x += 1
    end
    s - 1
end

@show projecteuler30p();
@btime projecteuler30p();

# %%
function projecteuler30()
    s = 0
    x = 0
    @nloops 6 i d -> 0:9 begin
        y = 0
        @nexprs 6 d -> y += (z = i_d^2; z^2*i_d)
        s += (x == y) * x
        x += 1
    end
    s - 1
end

@show projecteuler30();
@btime projecteuler30();

# %%
function projecteuler30x()
    s = 0
    x = 0
    @nloops 6 i d -> 0:9 d -> begin
        z = i_d^2
        y_d = z^2 * i_d
    end begin
        y = 0
        @nexprs 6 d -> y += y_d
        s += (x == y) * x
        x += 1
    end
    s - 1
end

@show projecteuler30x();
@btime projecteuler30x();

# %%
[(k, k * 9^4, log10(k * 9^4)+1) for k in 1:9]

# %%
[(k, k * 9^5, log10(k * 9^5)+1) for k in 1:9]

# %%
[(k, k * 9^10, log10(k * 9^10)+1) for k in 1:12]

# %%
[(m, findfirst(k -> k > log10(k*10.0^m)+1, 1:2m)-1) for m in 2:100]

# %%
?@nloops

# %%
(@macroexpand @nloops 3 i d -> 0:9 d -> y_d = i_d^2 begin
    y = 0
    @nexprs 3 d -> y += y_d
end)|> Base.remove_linenums!

# %%
@ntuple 10 i -> (i - 1) ^ 3

# %%
(@macroexpand @ntuple 10 i -> (i - 1) ^ 3) |> Base.remove_linenums!

# %%
