# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# $$
# R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}R = \frac{8\pi G}{c^4}T_{\mu\nu}
# $$

# %%
# https://nbviewer.jupyter.org/urls/gist.githubusercontent.com/terasakisatoshi/e4af4f618161d87cdc288285585527c5/raw/3a50de35748be12ab6461156440685052aadfcd2/hexagon_something.ipynb

using HTTP, JSON
function display_tweet(link)
    api = "https://publish.twitter.com/oembed?url=$link"
    r = response = HTTP.request("GET", api);
    j = JSON.parse(String(r.body))
    HTML(j["html"])
end
display_tweet("https://twitter.com/fujitapiroc1964/status/1436641682234114055")

# %%
