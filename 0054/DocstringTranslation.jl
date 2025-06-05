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
#     display_name: Julia current stable release
#     language: julia
#     name: julia
# ---

# %% colab={"base_uri": "https://localhost:8080/"} id="ftyR8UfgFA1y" outputId="386de2cc-9ca7-46a0-a6ba-95d24f68f87a"
# Google Colabと自分のパソコンの両方で使えるようにするための工夫

using Pkg

"""すでにPkg.add済みのパッケージのリスト (高速化のために用意)"""
_packages_added = [sort!(readdir(Sys.STDLIB));
    [info.name for (uuid, info) in Pkg.dependencies() if info.is_direct_dep]]

"""_packages_added内にないパッケージをPkg.addする"""
add_pkg_if_not_added_yet(pkg) = if !(pkg in _packages_added)
    println(stderr, "# $(pkg).jl is not added yet, so let's add it.")
    Pkg.add(pkg)
end

"""expr::Exprからusing内の`.`を含まないモジュール名を抽出"""
function find_using_pkgs(expr::Expr)
    pkgs = String[]
    function traverse(expr::Expr)
        if expr.head == :using
            for arg in expr.args
                if arg.head == :. && length(arg.args) == 1
                    push!(pkgs, string(arg.args[1]))
                elseif arg.head == :(:) && length(arg.args[1].args) == 1
                    push!(pkgs, string(arg.args[1].args[1]))
                end
            end
        else
            for arg in expr.args arg isa Expr && traverse(arg) end
        end
    end
    traverse(expr)
    pkgs
end

"""必要そうなPkg.addを追加するマクロ"""
macro autoadd(expr)
    pkgs = find_using_pkgs(expr)
    :(add_pkg_if_not_added_yet.($(pkgs)); $expr)
end

# %% colab={"base_uri": "https://localhost:8080/"} id="ztGw3D2i2Gaa" outputId="1b03c348-c7f5-4cab-c0f0-8a2063608625"
macro usingdocstringtranslation()
    quote
        translationurl = "https://github.com/AtelierArith/DocStrBankExperimental.jl/releases/download/full-ja/translation.zip"
        #translationurl = "https://github.com/AtelierArith/DocStrBankExperimental.jl/releases/download/stats-ja/translation.zip"
        scratchspaces_dir = joinpath(DEPOT_PATH[1], "scratchspaces", "d404e13b-1f8e-41a5-a26a-0b758a0c6c97")
        translation_dir = joinpath(scratchspaces_dir, "translation")
        if isdir(translation_dir)
            println(stderr, "Directory $translation_dir already exists.")
            println(stderr, "To download and extract translation.zip again, the directory must be deleted.")
        else
            println(stderr, "translation.zip to be downloaded will be extracted into $translation_dir.")
            run(`wget --no-verbose $(translationurl) -P $(scratchspaces_dir)`)
            run(`unzip -q $(joinpath(scratchspaces_dir, "translation.zip")) -d $(scratchspaces_dir)`)
        end
        @autoadd using DocstringTranslation
    end
end

# Google Colabで次の行の実行には3分秒程度かかる。
@usingdocstringtranslation; @switchlang! :ja

# %% colab={"base_uri": "https://localhost:8080/", "height": 1000} id="dAT1UTlgFxQV" outputId="2a9effed-e96c-4c28-a11b-45a94395fc77"
?range

# %% colab={"base_uri": "https://localhost:8080/", "height": 1000} id="pv_opZbtF27t" outputId="fd88bd7e-59a7-499c-accb-94b9ff645d7f"
??range

# %% colab={"base_uri": "https://localhost:8080/"} id="-CvHtpJwJ56y" outputId="e718c789-cedb-46e3-df7c-e5ab80b8c64d"
range(length=10)

# %% colab={"base_uri": "https://localhost:8080/"} id="1F3kNq2uSmu5" outputId="a20bf429-8087-4777-fba4-cf1f9d251503"
"""2倍する関数"""
double(x) = 2x

# %% colab={"base_uri": "https://localhost:8080/", "height": 81} id="WnwYsHBmD4in" outputId="e3ad75b3-6386-4552-a9be-887911d13c2f"
?double

# %% colab={"base_uri": "https://localhost:8080/"} id="q2Uv5o2kFLMD" outputId="3d7f5aa5-7393-48e9-e0f4-5b7eed10acc1"
@autoadd begin
using LinearAlgebra
using Distributions
using StatsPlots
using StatsBase
using StatsFuns
end

# %% colab={"base_uri": "https://localhost:8080/", "height": 792} id="W1-t8guzpUgG" outputId="375e2909-3255-4361-f93d-12629388125b"
?dot

# %% colab={"base_uri": "https://localhost:8080/", "height": 306} id="AjUckfbrHwtu" outputId="a3ae93da-22e5-42d3-db6d-8a6640e415b3"
?Gamma

# %% colab={"base_uri": "https://localhost:8080/", "height": 517} id="V3dIYTyJFbNF" outputId="452120fa-b62d-4ec3-e2df-79fa33efece2"
?density

# %% colab={"base_uri": "https://localhost:8080/", "height": 165} id="ALDiZm3RFc0S" outputId="e31db258-ab32-47d7-b908-f0e12c197cab"
?ecdf

# %% colab={"base_uri": "https://localhost:8080/", "height": 200} id="NKwp42CYFeQJ" outputId="0a001c5b-1f46-4a20-9ae8-adee4b80a338"
?logit

# %% id="zUVIRUhgFgDj"
