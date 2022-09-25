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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %%
module Humans

Base.@kwdef mutable struct Human
    name::String = "名無し"
    age::Int = 42
end

name(h::Human) = h.name
age(h::Human) = h.age
setname!(h::Human, name) = setfield!(h, :name, name)
setage!(h::Human, age) = setfield!(h, :age, age)

function hello(h::Human)
    "こんにちは。$(name(h))です。年齢は$(age(h))歳です。"
end

end

# %%
methods(Humans.Human)

# %%
nanashi = Humans.Human()

# %%
Humans.hello(nanashi)

# %%
Humans.setname!(nanashi, Humans.name(nanashi) * "の権平")
Humans.setage!(nanashi, Humans.age(nanashi) + 100)
Humans.hello(nanashi)

# %%
julia = Humans.Human(name="Julia")

# %%
morty = Humans.Human("Morty", 12)

# %%
ss = SubString("12しんじ34", 3, 9)

# %%
typeof(ss)

# %%
shinji = Humans.Human(name=ss, age=14.0)

# %%
typeof(Humans.name(shinji)), typeof(Humans.age(shinji))

# %%
module Humans

Base.@kwdef mutable struct Human{Tname, Tage}
    name::Tname = "名無し"
    age::Tage = 42
end

name(h::Human) = h.name
age(h::Human) = h.age
setname!(h::Human, name) = setfield!(h, :name, name)
setage!(h::Human, age) = setfield!(h, :age, age)

function hello(h::Human)
    "こんにちは。$(name(h))です。年齢は$(age(h))際です。"
end

end

# %%
methods(Humans.Human)

# %%
methods(Humans.Human{String, Int})

# %%
nanashi = Humans.Human()

# %%
Humans.hello(nanashi)

# %%
Humans.setname!(nanashi, Humans.name(nanashi) * "の権平")
Humans.setage!(nanashi, Humans.age(nanashi) + 100)
Humans.hello(nanashi)

# %%
julia = Humans.Human(name="Julia")

# %%
morty = Humans.Human("Morty", 12)

# %%
ss = SubString("12しんじ34", 3, 9)

# %%
typeof(ss)

# %%
shinji = Humans.Human(name=ss, age=14.0)

# %%
typeof(Humans.name(shinji)), typeof(Humans.age(shinji))

# %%
