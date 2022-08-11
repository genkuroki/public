### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ c0a7fcb0-84af-11eb-193f-155794030423
using MonteCarloMeasurements, Plots, PlutoUI

# ╔═╡ d8fc78e0-84af-11eb-0ae1-d9af9ccaf654
function symplectic_euler(g, f, v0, u0, tspan, param, dt)
	u_prev, v_prev = promote(u0, v0)
	t_prev, tmax = tspan
	t, u, v = [t_prev], [u_prev], [v_prev]
	while t_prev < tmax
		v_next = v_prev + g(u_prev, v_prev, param)*dt
		u_next = u_prev + f(u_prev, v_next, param)*dt
		t_next = t_prev + dt
		push!(t, t_next)
		push!(u, u_next)
		push!(v, v_next)
		t_prev, u_prev, v_prev = t_next, u_next, v_next
	end
	t, u, v
end

# ╔═╡ b13930a0-84be-11eb-2555-1925380667dd
md"""
__Van der Pol equation:__

``
\displaystyle\quad
\frac{du}{dt} = v, \quad
\frac{dv}{dt} = -\omega^2 u + \mu(1-u^2)v, \quad
u(0) = u_0, \quad
v(0) = v_0
``
"""

# ╔═╡ 3c6b2f20-84b0-11eb-36ee-a315d054a9c3
begin
	ω_slider = @bind ω Slider(0.1:0.1:5; default=2.0, show_value=true)
	μ_slider = @bind μ Slider(-1:0.1:2; default=0.3, show_value=true)
	u0_slider = @bind u00 Slider(-2:0.1:2; default=1.0, show_value=true)
	v0_slider = @bind v00 Slider(-2:0.1:2; default=0.0, show_value=true)
	tmax_slider = @bind tmax Slider(1.0:100; default=30.0, show_value=true)
	dt_slider = @bind dt Slider(0.01:0.01:0.1; default=0.05, show_value=true)
	ε_slider = @bind ε Slider(0:0.02:0.5; default=0.3, show_value=true)
	
	md"""
	ω: $ω_slider　μ: $μ_slider
	
	u₀: $u0_slider　v₀: $v0_slider
	
	tₘₐₓ: $tmax_slider　Δt: $dt_slider
	
	ε: $ε_slider
	"""
end

# ╔═╡ 77bfc860-84b0-11eb-3280-2b63049757be
begin 
	param = (ω, μ)
	f(u, v, param) = ((ω, μ) = param; v)
	g(u, v, param) = ((ω, μ) = param; -ω^2*u + μ*(1 - u^2)*v)
	u0, v0, tspan = u00 ± ε, v00, (0.0, tmax)
	t, u, v = symplectic_euler(g, f, v0, u0, tspan, param, dt)
	UV = plot(mean.(u), mean.(v); label="", c=:black, xlabel="u", ylabel="v", lw=0.5)
	U = plot(t, u; label="", c=:blue, xlabel="t", ylabel="u")
	V = plot(t, v; label="", c=:red,  xlabel="t", ylabel="v")
	layout = @layout[a{0.4w} grid(2,1)]
	plot(UV, U, V; size=(680, 300), layout, guidefontsize=7, tickfontsize=5)
end

# ╔═╡ Cell order:
# ╟─c0a7fcb0-84af-11eb-193f-155794030423
# ╟─d8fc78e0-84af-11eb-0ae1-d9af9ccaf654
# ╟─b13930a0-84be-11eb-2555-1925380667dd
# ╟─3c6b2f20-84b0-11eb-36ee-a315d054a9c3
# ╟─77bfc860-84b0-11eb-3280-2b63049757be
