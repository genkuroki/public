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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
using Nemo

# %%
F4, a = finite_field(2, 2, "a")

# %%
F64, b = finite_field(2, 6, "b")

# %%
F64(a)

# %%
collect_F4 = collect(F4)

# %%
F4_embedded_in_F64 = F64.(collect_F4)

# %%
a^2

# %%
b^6

# %%
ENV["LIES"] = 1000
p = 2
[((k, x) = finite_field(p, e, "x"); (p, e, x^e)) for e in 1:100]

# %%
p = 2
[(e = 2^k; (k, x) = finite_field(p, e, "x"); (p, e, x^e)) for k in 1:11]

# %%
p = 2
[(e = 2^k-1; (k, x) = finite_field(p, e, "x"); (p, e, x^e)) for k in 1:12]

# %%
p = 2
[(e = 2^k+1; (k, x) = finite_field(p, e, "x"); (p, e, x^e)) for k in 1:12]

# %%
F, a = finite_field(2, 1, :a)

# %%
A, x = polynomial_ring(F, :x)

# %%
is_irreducible(x^2+x+1)

# %%
is_irreducible(x^4+x+1)

# %%
ENV["LINES"] = 1000
F, a = finite_field(2, 1, :a)
A, x = polynomial_ring(F, :x)
[(p=2, e=e) for e in 2:100 if is_irreducible(x^e + x + 1)]

# %%
F, a = finite_field(2, 1, :a)
A, x = polynomial_ring(F, :x)
[(p=2, k=k, e=2^k) for k in 1:16 if (e = 2^k; is_irreducible(x^e + x + 1))]

# %%
F, a = finite_field(2, 1, :a)
A, x = polynomial_ring(F, :x)
[(p=2, k=k, e=2^k-1) for k in 1:14 if (e = 2^k-1; is_irreducible(x^e + x + 1))]

# %%
F, a = finite_field(2, 1, :a)
A, x = polynomial_ring(F, :x)
[(p=2, k=k, e=2^k+1) for k in 1:14 if (e = 2^k+1; is_irreducible(x^e + x + 1))]

# %%
