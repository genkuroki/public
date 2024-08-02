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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# https://mathworld.wolfram.com/q-PolygammaFunction.html
#
# The $q$-digamma function $\psi_q(z)$ can be defined by
#
# $$
# \psi_q(z) = \log(1 - q) + \log(q) \sum_{n=0}^\infty \frac{1}{q^{-(n+z)} - 1}.
# $$
#
# Let us plot the graph of $\displaystyle \tilde{\psi}_q(z) = \sum_{n=0}^\infty \frac{1}{q^{-(n+z)} - 1}$.

# %%
using CUDA
using Plots
default(fmt=:png)
using Printf

approxlogabs(x; a=0.1) = (abs(x)^a - 1)/a

qdigamma(q, z; L=100) = sum(n -> 1/(q ^ (-(n+z)) - 1), 0:L)

function qdigammadot(q, z; L=100)
    w = zero(q)
    for n in 0:L
        w .+= 1 ./ (q .^ (-(n+z)) .- 1)
    end
    w
end

function plot_qdigamma(;
        z=1,
        L=100,
        f=approxlogabs, 
        size=(500, 520),
        color=:glasbey_bw_minc_20_hue_150_280_n256,
        q=complex.(range(-1.1, 1.1, 400)', range(-1.1, 1.1, 400)),
        title="z = $(@sprintf "%5.2f" z),  f = $f", 
        kwargs...
    )
    w = qdigammadot(q, z; L)
    a = f.(w)
    x, y = real(q[1,:]), imag(q[:,1])
    heatmap(Array(x), Array(y), Array(a); colorbar=false, color, size, title, kwargs...)
end

# %%
plot(log, 0.0001, 2; label="log")
plot!(approxlogabs, 0.0001, 2; label="approxlogabs")

# %%
q = complex.(range(-1.1, 1.1, 400)', range(-1.1, 1.1, 400))
@show typeof(q)
@time a = approxlogabs.(qdigammadot(q, 2.0; L=1000))
@time plot_qdigamma(; z=0.1, L=1000)

# %%
q = complex.(range(-1.1, 1.1, 400)', range(-1.1, 1.1, 400))
q_cu = CuArray(q)
@show typeof(q_cu)
@time a = approxlogabs.(qdigammadot(q_cu, 2.0; L=1000))
@time a = approxlogabs.(qdigammadot(q_cu, 2.0; L=1000))
@time a = approxlogabs.(qdigammadot(q_cu, 2.0; L=1000))
@time plot_qdigamma(; z=0.1, L=1000, q=q_cu)

# %%
q = q_cu; # Comment out this line if you don't use CUDA.

# %%
@time plot_qdigamma(; z=0.1, f=angle, L=1000, q)

# %%
@time plot_qdigamma(; z=1, q)

# %%
@time plot_qdigamma(; z=1, f=angle, q)

# %%
@time plot_qdigamma(; z=2, q)

# %%
@time plot_qdigamma(; z=2, f=angle, q)

# %%
@time anim = @animate for z in [fill(0.1, 10); 0.1:0.02:2; fill(2.0, 10)]
    plot_qdigamma(; z=z, clim=(-5.65, 9.95), q)
end
gif(anim, "qdigamma.gif")

# %%
@time anim = @animate for z in [fill(0.1, 10); 0.1:0.02:2; fill(2.0, 10)]
    plot_qdigamma(; z=z, f=angle, q)
end
gif(anim, "qdigamma2.gif")

# %%
