---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.3
  kernelspec:
    display_name: Julia 1.7.0
    language: julia
    name: julia-1.7
---

```julia
using Distributions
using Plots

using Base64
showimg(mime, fn; tag="img") = open(fn) do f
    base64 = base64encode(f)
    display("text/html", """<$tag src="data:$mime;base64,$base64" />""")
end

plot(sin; size=(200, 150));
```

```julia
dist1 = MixtureModel([Normal(0, 0.1), Normal(1, 0.1)], [2/3, 1/3])
```

```julia
dist2 = MixtureModel([Normal(0, 0.05), Normal(1, 0.05)], [2/3, 1/3])
```

```julia
P = plot(; fontfamily="Computer Modern")
plot!(x -> pdf(dist1, x), -0.5, 1.5; label="\$\\mu_1=0,\\; \\mu_2=1,\\; \\sigma=0.1\$")
plot!(x -> pdf(dist2, x), -0.5, 1.5; label="\$\\mu_1=0,\\; \\mu_2=1,\\; \\sigma=0.05\$")
savefig(P, "test100.png")
showimg("image/png", "test100.png"; tag="img width=$(P.attr[:size][1])")
```

```julia
P.attr[:dpi] = 600
savefig(P, "test600.png")
showimg("image/png", "test600.png"; tag="img width=$(P.attr[:size][1])")
```

```julia

```
