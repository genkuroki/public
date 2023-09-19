# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# [左右対称で平均と分散と尖度が一致するがMann-WhitneyのU検定を使ってはいけないシンプルな具体例.ipynb](https://github.com/genkuroki/public/blob/main/0041/%E5%B7%A6%E5%8F%B3%E5%AF%BE%E7%A7%B0%E3%81%A7%E5%B9%B3%E5%9D%87%E3%81%A8%E5%88%86%E6%95%A3%E3%81%A8%E5%B0%96%E5%BA%A6%E3%81%8C%E4%B8%80%E8%87%B4%E3%81%99%E3%82%8B%E3%81%8CMann-Whitney%E3%81%AEU%E6%A4%9C%E5%AE%9A%E3%82%92%E4%BD%BF%E3%81%A3%E3%81%A6%E3%81%AF%E3%81%84%E3%81%91%E3%81%AA%E3%81%84%E3%82%B7%E3%83%B3%E3%83%97%E3%83%AB%E3%81%AA%E5%85%B7%E4%BD%93%E4%BE%8B.ipynb) に名前を変えた→[nbviewer](https://nbviewer.org/github/genkuroki/public/blob/main/0041/%E5%B7%A6%E5%8F%B3%E5%AF%BE%E7%A7%B0%E3%81%A7%E5%B9%B3%E5%9D%87%E3%81%A8%E5%88%86%E6%95%A3%E3%81%A8%E5%B0%96%E5%BA%A6%E3%81%8C%E4%B8%80%E8%87%B4%E3%81%99%E3%82%8B%E3%81%8CMann-Whitney%E3%81%AEU%E6%A4%9C%E5%AE%9A%E3%82%92%E4%BD%BF%E3%81%A3%E3%81%A6%E3%81%AF%E3%81%84%E3%81%91%E3%81%AA%E3%81%84%E3%82%B7%E3%83%B3%E3%83%97%E3%83%AB%E3%81%AA%E5%85%B7%E4%BD%93%E4%BE%8B.ipynb).
