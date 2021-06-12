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

<!-- #region -->
See https://github.com/ranaroussi/yfinance

To install, run

```bash
pip install yfinance
```
<!-- #endregion -->

```julia
using PyCall
yf = pyimport("yfinance")
```

```julia
msft = yf.Ticker("MSFT")
```

```julia
msft.info
```

```julia
hist = msft.history(period="max")
```

```julia
msft.actions
```

```julia

```
