# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:hydrogen
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import sys
sys.version

# %% [markdown]
# https://zenn.dev/turing_motors/articles/8b9a2c4d3e8882

# %%
import time

xmin = -1.75
xmax = 0.75
width = 4096
ymin = -1.25
ymax = 1.25
height = 4096
max_iter = 500

def mandelbrot_kernel(c):
    z = c
    for i in range(max_iter):
        z = z * z + c
	# zが閾値を超えたら終了します
        if abs(z) > 2:
            return i
    return max_iter

def compute_mandelbrot(image):
    dx = (xmax - xmin) / width
    dy = (ymax - ymin) / height
		
    # 各ピクセルごとに複素数を計算します
    for j in range(height):
        for i in range(width):
            y = ymin + j * dy
            x = xmin + i * dx
            image[j][i] = mandelbrot_kernel(complex(x, y))
    return image

image = [[0 for _ in range(width)] for _ in range(height)]
t0 = time.time()
image = compute_mandelbrot(image)
t1 = time.time()
print(f"elapsed time: {(t1-t0) * 1000:.3f} ms")

# %%
