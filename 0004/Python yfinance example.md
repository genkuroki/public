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
    display_name: Python 3
    language: python
    name: python3
---

<!-- #region -->
See https://github.com/ranaroussi/yfinance

To install, run

```bash
pip install yfinance
```
<!-- #endregion -->

```python
import yfinance as yf
```

```python
msft = yf.Ticker("MSFT")
```

```python
msft.info
```

```python
hist = msft.history(period="max")
```

```python
msft.actions
```

```python

```
