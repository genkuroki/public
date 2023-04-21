% -*- coding: utf-8 -*-
% ---
% jupyter:
%   jupytext:
%     formats: ipynb,m:hydrogen
%     text_representation:
%       extension: .m
%       format_name: hydrogen
%       format_version: '1.3'
%       jupytext_version: 1.10.3
%   kernelspec:
%     display_name: Wolfram Language 13.1
%     language: Wolfram Language
%     name: wolframlanguage13.1
% ---

% %% tags=[]
RegionPlot[ImplicitRegion[x y > x^4 + y^4, {x, y}]]

% %%
NIntegrate[1, {x, y} ∈ ImplicitRegion[x y > x^4 + y^4, {x, y}]]/2

% %%
N[Beta[1/2, 1/2]/8]

% %%
NIntegrate[1, {x, y} ∈ ImplicitRegion[x y > x^6 + y^6, {x, y}]]/2

% %%
N[Beta[1/4, 1/4]/12]

% %%
NIntegrate[1, {x, y} ∈ ImplicitRegion[x y > x^8 + y^8, {x, y}]]/2

% %%
N[Beta[1/6, 1/6]/16]

% %%
