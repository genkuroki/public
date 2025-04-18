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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
import Pkg
packages_added = [info.name for (uuid, info) in Pkg.dependencies() if info.is_direct_dep]
for pkg in ("Distributions", "StatsFuns")
    if pkg in packages_added
        println("$(pkg).jl is already added.")
    else
        println("$(pkg).jl is not added yet, do let's add it.")
        Pkg.add(pkg)
    end
end

# %% colab={"base_uri": "https://localhost:8080/"} executionInfo={"elapsed": 3661, "status": "ok", "timestamp": 1744750950679, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="ZLDGc0s7fPyS" outputId="27b0d711-8953-4133-ff28-769e574889aa"
using Random
using Distributions
using StatsFuns
using Plots
default(fmt=:png, legendfontsize=12)
safediv(x, y) = x == 0 ? zero(x/y) : x/y

varstab_bin(p) = asin(2p - 1)

function expectval(f, bin::Binomial)
    sum(f(x) * pdf(bin, x) for x in support(bin))
end

function variance(f, bin::Binomial)
    mu = expectval(f, bin)
    sigma2 = expectval(x -> (f(x) - mu)^2, bin)
end

function variance_varstab_bin(n, p)
    variance(k -> sqrt(n) * varstab_bin(k/n), Binomial(n, p))
end

function variance_logit(n, p)
    variance(k -> 0 < k < n ? sqrt(n)/2 * logit(k/n) : 0.0, Binomial(n, p))
end

function pvalue_binomial_wald(k, n, p)
    phat = k/n
    z = safediv(phat - p, sqrt(phat*(1 - phat)/n))
    2ccdf(Normal(), abs(z))
end

function pvalue_binomial_score(k, n, p)
    phat = k/n
    z = safediv(phat - p, sqrt(p*(1 - p)/n))
    2ccdf(Normal(), abs(z))
end

function pvalue_binomial_varstab(k, n, p)
    phat = k/n
    correction = (2p-1)/(4n*sqrt(p*(1-p)))
    z = sqrt(n) * (varstab_bin(phat) - varstab_bin(p) - correction)
    2ccdf(Normal(), abs(z))
end

function alpha_erroer_rate(pvaluefunc, bin::Binomial; alpha=0.05)
    n, p = params(bin)
    f(k) = (pvaluefunc(k, n, p) < alpha)
    expectval(f, bin)
end

function plot_comparison_varstab_logit(; n = 10, f=trues(2), kwargs...)
    plot()
    f[1] && plot!(p -> variance_varstab_bin(n, p), 0, 1; label="var.stab.")
    f[2] && plot!(p -> variance_logit(n, p), 0, 1; label="logit")
    hline!([1]; c=:black, ls=:dot, label="")
    plot!(legend=:top)
    plot!(; kwargs...)
end

function plot_pvalue_functions(; k=1, n=10, p=0.0, kwargs...)
    @show k, n, p
    @show pvalue_binomial_varstab.(k, n, p)
    @show pvalue_binomial_score.(k, n, p)
    @show pvalue_binomial_wald.(k, n, p)
    plot()
    plot!(p -> pvalue_binomial_varstab(k, n, p), 0, 1; label="var.stab.")
    plot!(p -> pvalue_binomial_score(k, n, p), 0, 1; label="score", ls=:dash)
    plot!(p -> pvalue_binomial_wald(k, n, p), 0, 1; label="Wald", ls=:dot)
    plot!(xtick=0:0.1:1, ytick=0:0.05:1)
    title!("P-value functions for k=$k, n=$n")
    plot!(; kwargs...)
end

function plot_alpha_error_rates(; alpha=0.5, n=10, ps=range(0, 1, 1001), f=trues(3), kwargs...)
    plot()
    f[1] && plot!(ps, p -> alpha_erroer_rate(pvalue_binomial_varstab, Binomial(n, p); alpha); label="var.stab.", c=1)
    f[2] && plot!(ps, p -> alpha_erroer_rate(pvalue_binomial_score, Binomial(n, p); alpha); label="score", c=2)
    f[3] && plot!(ps, p -> alpha_erroer_rate(pvalue_binomial_wald, Binomial(n, p); alpha); label="wald", ls=:dot, c=3)
    hline!([alpha]; label="", ls=:dot, c=:black)
    plot!(xtick=0:0.1:1, ytick=0:0.05:1)
    plot!(legend=:top)
    title!("alpha error rates for n=$n, alpha=$alpha")
    plot!(; kwargs...)
end

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 4392, "status": "ok", "timestamp": 1744750955073, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="kxPMhMSwazYV" outputId="f8bd8b14-63b7-400f-aa95-75ac62a43a98"
plot_comparison_varstab_logit(; n = 10)

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 141, "status": "ok", "timestamp": 1744750955216, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="uiVXeniebX1N" outputId="7607ce09-ea3e-46cc-e28b-d85c232c85a4"
plot_comparison_varstab_logit(; n = 100)

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 85, "status": "ok", "timestamp": 1744750955301, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="qnRut2btb11C" outputId="5841c5cb-28e0-42de-b569-f3dec2f0d280"
plot_comparison_varstab_logit(; n = 100, ylim=(-0.1, 3.1))

# %% colab={"base_uri": "https://localhost:8080/", "height": 496} executionInfo={"elapsed": 827, "status": "ok", "timestamp": 1744750956130, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="_f0R_VLPcrFp" outputId="56e7829d-59c4-4b82-bf8d-1711da0a4467"
plot_pvalue_functions(; k=1, n=10, p=(0.0, 0.4))

# %% colab={"base_uri": "https://localhost:8080/", "height": 496} executionInfo={"elapsed": 138, "status": "ok", "timestamp": 1744750956269, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="ufAfnwDOc2J2" outputId="7aadfbc3-01c0-4e2b-86c2-5f7e29fe4ce9"
plot_pvalue_functions(; k=2, n=10, p=(0.0, 0.5))

# %% colab={"base_uri": "https://localhost:8080/", "height": 496} executionInfo={"elapsed": 319, "status": "ok", "timestamp": 1744750956589, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="Iawf4xzqdA13" outputId="741e3d45-1e69-429a-a9f3-f353891962e4"
plot_pvalue_functions(; k=20, n=100, p=0.285, xlim=(0.05, 0.4), xtick=0:0.05:1)

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 1414, "status": "ok", "timestamp": 1744750958006, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="BSzZu1Lbdz8B" outputId="32ce03a7-384b-4b78-f2e5-80c1648b9bbd"
plot_alpha_error_rates(; alpha=0.05, n=10)

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 177, "status": "ok", "timestamp": 1744750958184, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="mP7Yq2eiM4Vp" outputId="1e4b0de1-89a6-45cb-f1b2-e4f25d32c5d5"
plot_alpha_error_rates(; alpha=0.05, n=10, f=Bool[1,1,0], ylim=(-0.01, 0.35))

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 276, "status": "ok", "timestamp": 1744750958461, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="m1cU8g_xd-FS" outputId="cf1d7944-61fc-4e77-f848-3294067088b9"
plot_alpha_error_rates(; alpha=0.05, n=30)

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 201, "status": "ok", "timestamp": 1744750958663, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="zXhIIj6WNF8k" outputId="a65f9508-22bc-4462-8464-d5d7b71840e1"
plot_alpha_error_rates(; alpha=0.05, n=30, f=Bool[1,1,0], ylim=(-0.01, 0.35))

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 299, "status": "ok", "timestamp": 1744750958963, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="7PVsEPACqOMP" outputId="760f27de-4f28-4743-80eb-1e3e46344e83"
plot_alpha_error_rates(; alpha=0.05, n=100)

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 185, "status": "ok", "timestamp": 1744750959149, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="JPuhsfBEeXyp" outputId="dfdc1141-9762-4116-9573-fe591d571aab"
plot_alpha_error_rates(; alpha=0.05, n=100, f=Bool[1,1,0], ylim=(-0.01, 0.35))

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 332, "status": "ok", "timestamp": 1744750959483, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="DIVF2j6rNSJE" outputId="ae64b55c-3410-4490-b230-d6cb64f735df"
plot_alpha_error_rates(; alpha=0.05, n=100, f=Bool[1,1,0], ylim=(-0.003, 0.113), ytick=0:0.01:1)

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 423, "status": "ok", "timestamp": 1744750959907, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="yXWzoXfYejLz" outputId="9a2b18a2-b8e0-4b60-850a-4ec3206cb07b"
plot_alpha_error_rates(; alpha=0.05, n=300, f=Bool[1,1,0], ytick=0:0.01:1, ylim=(-0.003, 0.113))

# %% colab={"base_uri": "https://localhost:8080/", "height": 422} executionInfo={"elapsed": 300, "status": "ok", "timestamp": 1744750960209, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="gIj2JG20Ot5E" outputId="437b32cd-a56d-4822-9d80-ab17691de11a"
plot_alpha_error_rates(; alpha=0.05, n=300, f=Bool[0,1,0], ytick=0:0.01:1, legend=:bottom)

# %% executionInfo={"elapsed": 1, "status": "ok", "timestamp": 1744750960211, "user": {"displayName": "Gen Kuroki", "userId": "09579187526869050279"}, "user_tz": -540} id="MRi_r9b0gQUv"
