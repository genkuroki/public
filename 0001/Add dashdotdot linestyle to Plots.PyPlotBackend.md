---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Julia 1.7.0-DEV
    language: julia
    name: julia-1.7
---

```julia
using Plots
pyplot(fmt = :svg)

function testplot()
    P = plot(legend=:outerright)
    for i in 1:20
        plot!(fill(-i, 10); label="$i", ls=:auto, ylim=(-21, 0))
    end
    P
end

testplot()
```

```julia
# Add :dashdotdot linestyle to Plots.PyPlotBackend()

@eval Plots begin
    push!(_pyplot_style, :dashdotdot)

    const dic_py_linestyle = Dict(
        :none => " ",
        :solid => "-",
        :dash => "--",
        :dot => ":",
        :dashdot => "-.",
        :dashdotdot => (0, (6, 1.2, 1.5, 1.2, 1.5, 1.2)),
    )

    function py_linestyle(seriestype::Symbol, linestyle::Symbol)
        ls = get(dic_py_linestyle, linestyle, nothing)
        !isnothing(ls) && return ls
        @warn("Unknon linestyle $linestyle")
        return "-"
    end
end
```

```julia
testplot()
```

```julia
plot()
for (k, linestyle) in enumerate(Plots._pyplot_style[2:end])
    plot!(x -> sin(x) - k, -5, 5; label="sin(x) - $k", ylim=(-6.5, 1.0), linestyle)
end
plot!()
```

```julia

```
