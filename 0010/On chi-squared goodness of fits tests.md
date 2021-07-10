The Chi-square test can only be used if the distribution of `testValue`s is approximated by the Chi-square distribution.

In the standard goodness-of-fit test, the degree of freedom of the chi-square distribution should be equal to (the number of elements in the support of the discrete distribution to be tested) - 1.

However, the real support of `Binomial(100, 0.7)` is at most between 40 and 100. So, in this case, we can expect that the degree of freedom should be roughly around 60.

Let's check this in the following.

The mathematically idealized support of `DiscreteUniform(0, 100)` and `Binomial(100, 0.7)` are both 0:100, but the real support of the latter is much smaller than the former.

```julia
using Distributions, StatsPlots
P = plot(DiscreteUniform(0, 100); label="DiscreteUniform(0, 100)", xtick=0:10:100, ylim=(0, 0.013))
Q = plot(Binomial(100, 0.7); label="Binomial(100, 0.7)", legend=:topleft, xtick=0:10:100)
plot(P, Q; size=(800, 250))
```

![images/2021-07-10 (0).png](https://raw.githubusercontent.com/genkuroki/public/main/0010/images/2021-07-10%20(0).png)

```julia
using Distributions, Random, StatsPlots

# The test definitions.....
function computeDensity(data, supp)
    counts =  [count(i -> i==s,data) for s in supp]
    if length(data) > sum(counts)
        error("There are some data not in the support !")
    end
    return counts
end

"""
Modified version of `goodnessOfFitDiscrete` function in
https://discourse.julialang.org/t/whats-wrong-with-my-chi-squared-goodness-of-fits-tests/64334

* The old argument `support` is deleted.
* The new local variable `supp` is defined to be `support(f₀)`
"""
function goodnessOfFitDiscrete(data, f₀; compressedData=true, α=0.05, d=0)
    supp       = support(f₀)
    if !compressedData
        data   = computeDensity(data, supp)
    end
    K          = length(supp)
    N          = sum(data)
    p̂          = data ./ N
    df         = K - d - 1
    p0         = pdf.(f₀,supp)
    T          = N * sum((p̂[k] - p0[k])^2/p0[k] for k in 1:K)
    χDist      = Chisq(df)
    rejectedH₀ = T > quantile(χDist, 1-α)
    p_value    = 1 - cdf(χDist, T)
    return (testValue=T, threashold=quantile(χDist,1-α), rejectedH₀=rejectedH₀, p_value=p_value)
end

function repeat_tests(f₀; datasize = 10000, repetitions = 10000, α = 0.05, d = 0)
    testValue = zeros(repetitions)
    rejectedH₀ = falses(repetitions)
    data = rand(f₀, datasize)
    for rep in 1:repetitions
        rand!(f₀, data)
        out = goodnessOfFitDiscrete(data, f₀; compressedData=false, α, d)
        testValue[rep] = out.testValue
        rejectedH₀[rep] = out.rejectedH₀
    end
    α_real = sum(rejectedH₀)/repetitions
    (; α_real, testValue, rejectedH₀)
end

function plot_testValue_dist(f₀; datasize = 10000, repetitions = 10000, α = 0.05, d=0)
    α_real, testValue, rejectedH₀ = repeat_tests(f₀; datasize, repetitions, α, d)

    title = "$f₀   (real α = $α_real / nominal α = $α)"
    plot(; title, titlefontsize=10)

    xlim=(0, 2length(support(f₀)))
    a, b = xlim
    histogram!(testValue[a .≤ testValue .≤ b]; norm=true, alpha=0.3, label="testValue", xlim)

    df = length(support(f₀)) - d - 1
    plot!(Chisq(df), a, b; label="Chisq(df = $df)", ls=:dash, lw=2)
end
```

In the `DiscreteUniform(0, 100)` case, the distribution of `testValue`s is almost exactly equal to the chi-square distribution.

```julia
plot_testValue_dist(DiscreteUniform(0, 100))
```

![images/2021-07-10 (1).png](https://raw.githubusercontent.com/genkuroki/public/main/0010/images/2021-07-10%20(1).png)

But, in the `Binomial(100, 0.7)` case, this is not the case.  The real (or empirical) α value is much smaller than the nominal α value.

```julia
plot_testValue_dist(Binomial(100, 0.7))
```

![images/2021-07-10 (2).png](https://raw.githubusercontent.com/genkuroki/public/main/0010/images/2021-07-10%20(2).png)

When the degree of freedom is set to 60, the empirical α value is close to the nominal α value 0.05.

```julia
plot_testValue_dist(Binomial(100, 0.7); d = 40)
```

![images/2021-07-10 (3).png](https://raw.githubusercontent.com/genkuroki/public/main/0010/images/2021-07-10%20(3).png)
