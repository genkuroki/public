# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
module LifeGame

using Plots
using ProgressMeter

struct Board{T}
    state::T
    tmp::T
end
Board(state) = Board(state, similar(state))
function fivexfive(n=200)
    state = zeros(Int8, n, n)
    m = n ÷ 2
    state[m-2:m+2, m-2:m+2] .= [1 1 1 0 1; 1 0 0 0 0; 0 0 0 1 1; 0 1 1 0 1; 1 0 1 0 1]
    Board(state)
end
randboard(n=200) = Board(rand(Int8[0, 1], n, n))

P(i, m) = ifelse(i == m, 1, i+1)
Q(i, m) = ifelse(i == 1, m, i-1)
function _update!(v, u)
    m, n = size(u)
    @inbounds for j in 1:n, i in 1:m
        i₊, i₋, j₊, j₋ = P(i, m), Q(i, m), P(j, n), Q(j, n)
        N = u[i₋,j₋] + u[i,j₋] + u[i₊,j₋] + u[i₋,j] + u[i₊,j] + u[i₋,j₊] + u[i,j₊] + u[i₊,j₊]
        v[i,j] = N == 3 ? 1 : (N == 2 && !iszero(u[i,j])) ? 1 : 0
    end
end

function update!(board::Board, niters=1)
    state, tmp = board.state, board.tmp
    for _ in 1:niters÷2
        _update!(tmp, state)
        _update!(state, tmp)
    end
    if isodd(niters)
        _update!(tmp, state)
        state .= tmp
    end
end

function gif!(board::Board, niters, nskips=1; gifname="life.gif", fps=20, size=(240, 240))
    prog = Progress(niters, 0)
    anim = @animate for t in 1:niters
        heatmap(board.state; size, colorbar=false, ticks=false, axis=false, frame=false)
        update!(board, nskips)
        next!(prog)
    end
    gif(anim, gifname; fps)
end

end

# %%
boardrandom = LifeGame.randboard()
LifeGame.gif!(boardrandom, 500; gifname="liferandom.gif")

# %%
board5x5 = LifeGame.fivexfive()
LifeGame.gif!(board5x5, 2000; gifname="life5x5.gif", fps=60)

# %%
state = Int8[
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 1 1 0 0 0 0
    1 0 1 0 0 0 0
    0 0 1 0 0 0 0
]
board = LifeGame.Board(state)
board.state

# %%
LifeGame.update!(board, 12)
board.state

# %%
LifeGame.gif!(board, 27; gifname="lifeglider.gif", fps=5)

# %%
