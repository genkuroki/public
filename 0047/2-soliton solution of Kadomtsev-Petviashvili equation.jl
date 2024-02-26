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
#     display_name: Julia 1.10.1
#     language: julia
#     name: julia-1.10
# ---

# %%
using Plots
default(fmt=:png)

using ForwardDiff
using SymPy

# %%
sol1(x, y, t; u=1, v=0, c=0) = (u-v)^2/2 * 1/cosh((c + (u-v)*x + (u^2-v^2)*y + (u^3-v^3)*t)/2)^2

xs = range(-10, 10, 100)
ys = range(-10, 10, 100)

@gif for t in range(-12, 12, 100)
    surface(xs, ys, (x,y)->sol1(x,y,t; v=-0.5); camera=(30, 75), colorbar=false)
    plot!(size=(500, 600))
    plot!(zlim=(0, 1.3))
end

# %%
xs = range(-5, 5, 100)
ys = range(-5, 5, 100)

@gif for t in range(-7, 7, 100)
    surface(xs, ys, (x,y)->sol1(x,y,t; v=-1); camera=(30, 75), colorbar=false)
    plot!(size=(500, 600))
    plot!(zlim=(0, 2))
end

# %%
xs = range(-5, 5, 100)

@gif for t in range(-5, 5, 100)
    plot(xs, x->sol1(x,0,t; v=-1); label="")
end

# %%
tau1(x, y, t; u=1, v=0, c=0) = 1 + exp(c + (u-v)*x + (u^2-v^2)*y + (u^3-v^3)*t)

@syms x y t u v c

tau1(x, y, t; u, v, c) |>display

expr1 = 2SymPy.diff(log(tau1(x, y, t; u, v, c)), x, x) |> simplify
expr1 |> display

sol1(x, y, t; u, v, c)

# %%
tau2(x, y, t; u1=1, v1=0, c1=0, u2=-1, v2=0, c2=0) = 1 + exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) + exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t) + (u1-u2)*(v1-v2)/((u1-v2)*(v1-u2)) * exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) * exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t)

@syms x y t

expr2 = 2SymPy.diff(log(tau2(x, y, t)), x, x) |> simplify
expr2 |> display

code2 = string(expr2)
code2 |> display

Meta.parse("sol2(x, y, t) = " * code2) |> eval

xs = range(-10, 10, 100)
ys = range(-10, 10, 100)

@gif for t in range(-10, 10, 100)
    surface(xs, ys, (x,y)->sol2(x,y,t); camera=(30, 75), colorbar=false)
    plot!(size=(500, 600))
end

# %%
tau2(x, y, t; u1=1, v1=0, c1=0, u2=-0.7, v2=0, c2=0) = 1 + exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) + exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t) + (u1-u2)*(v1-v2)/((u1-v2)*(v1-u2)) * exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) * exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t)

@syms x y t

expr2 = 2SymPy.diff(log(tau2(x, y, t)), x, x) |> simplify
expr2 |> display

code2 = string(expr2)
code2 |> display

Meta.parse("sol2(x, y, t) = " * code2) |> eval

xs = range(-10, 10, 100)
ys = range(-10, 10, 100)

@gif for t in range(-10, 10, 100)
    surface(xs, ys, (x,y)->sol2(x,y,t); camera=(30, 75), colorbar=false)
    plot!(size=(500, 600))
end

# %%
tau2(x, y, t; u1=1, v1=0, c1=0, u2=1, v2=2, c2=0) = 1 + exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) + exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t) + (u1-u2)*(v1-v2)/((u1-v2)*(v1-u2)) * exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) * exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t)

@syms x y t

expr2 = 2SymPy.diff(log(tau2(x, y, t)), x, x) |> simplify
expr2 |> display

code2 = string(expr2)
code2 |> display

Meta.parse("sol2(x, y, t) = " * code2) |> eval

xs = range(-20, 20, 200)
ys = range(-20, 20, 200)

@gif for t in range(-5, 5, 100)
    surface(xs, ys, (x,y)->sol2(x,y,t); camera=(30, 75), colorbar=false)
    plot!(size=(500, 600))
    plot!(zlim=(0, 2))
end

# %%
tau2(x, y, t; u1=1, v1=-1, c1=0, u2=1.5, v2=-1.5, c2=0) = 1 + exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) + exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t) + (u1-u2)*(v1-v2)/((u1-v2)*(v1-u2)) * exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) * exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t)

@syms x y t

expr2 = 2SymPy.diff(log(tau2(x, y, t)), x, x) |> simplify
expr2 |> display

code2 = string(expr2)
code2 |> display

Meta.parse("sol2(x, y, t) = " * code2) |> eval

xs = range(-10, 10, 200)
ys = range(-10, 10, 200)

@gif for t in range(-5, 5, 100)
    surface(xs, ys, (x,y)->sol2(x,y,t); camera=(30, 75), colorbar=false)
    plot!(size=(500, 600))
    plot!(zlim=(0, 4.5))
end

# %%
xs = range(-10, 10, 300)

@gif for t in range(-6, 6, 100)
    plot(xs, x->sol2(x,0,t); label="")
    plot!(ylim=(-0.1, 4.5))
end

# %%
tau2(x, y, t; u1=1, v1=-1, c1=0, u2=1.5, v2=-1.2, c2=0) = 1 + exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) + exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t) + (u1-u2)*(v1-v2)/((u1-v2)*(v1-u2)) * exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) * exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t)

@syms x y t

expr2 = 2SymPy.diff(log(tau2(x, y, t)), x, x) |> simplify
expr2 |> display

code2 = string(expr2)
code2 |> display

Meta.parse("sol2(x, y, t) = " * code2) |> eval

xs = range(-20, 20, 200)
ys = range(-10, 10, 200)

@gif for t in range(-10, 10, 100)
    surface(xs, ys, (x,y)->sol2(x,y,t); camera=(30, 75), colorbar=false)
    plot!(size=(500, 600))
    plot!(zlim=(0, 4.5))
end

# %% tags=[]
tau2(x, y, t; u1=1, v1=-1, c1=0, u2=1.2, v2=-1.5, c2=0) = 1 + exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) + exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t) + (u1-u2)*(v1-v2)/((u1-v2)*(v1-u2)) * exp(c1 + (u1-v1)*x + (u1^2-v1^2)*y + (u1^3-v1^3)*t) * exp(c2 + (u2-v2)*x + (u2^2-v2^2)*y + (u2^3-v2^3)*t)

@syms x y t

expr2 = 2SymPy.diff(log(tau2(x, y, t)), x, x) #|> simplify
#expr2 |> display

code2 = string(expr2)
code2 |> display

Meta.parse("sol2(x, y, t) = " * code2) |> eval

xs = range(-20, 20, 200)
ys = range(-10, 10, 200)

@gif for t in range(-10, 10, 100)
    surface(xs, ys, (x,y)->sol2(x,y,t); camera=(30, 75), colorbar=false)
    plot!(size=(500, 600))
    plot!(zlim=(0, 4.5))
end

# %%
