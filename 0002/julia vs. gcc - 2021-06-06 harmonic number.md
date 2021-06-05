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
versioninfo()
```

```julia
versioninfo()
println()

function f(x, T=Float64)
    n = 1
    s = one(T)
    while s < x
        n += 1
        s += inv(T(n))
    end
    n, s
end

@time f(21)
```

```julia
using BenchmarkHistograms
@benchmark f(21)
```

```julia
run(`gcc --version`)
flush(stdout)

C_code = raw"""
long long f(double x) {
    long long n = 1;
    double s = 1.0;
    while (s < x) {
        n++;
        s += 1.0 / (double) n;
    }
    return n;
}
"""
display("text/markdown", "```C\n"*C_code*"\n```")

using Libdl
libname = tempname()
libname_dl = libname * "." * Libdl.dlext

open(`gcc -Wall -O3 -march=native -xc -shared -o $libname_dl -`, "w") do f
     print(f, C_code)
end
run(`ls -l $libname_dl`)
println()

f_gcc(x::Float64) = @ccall libname.f(x::Float64)::Int64
@time f_gcc(21.0)
```

```julia
@benchmark 6000125006293f_gcc(21.0)
```

```julia
0.67/740461601*6000125006293/60^2
```

```julia

```
