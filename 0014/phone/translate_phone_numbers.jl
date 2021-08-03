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
# https://discourse.julialang.org/t/help-to-get-my-slow-julia-code-to-run-as-fast-as-rust-java-lisp/65741

# %%
# download dictionary
if !isfile("dictionary.txt")
    dictionary_url = "https://raw.githubusercontent.com/renatoathaydes/prechelt-phone-number-encoding/julia/dictionary.txt"
    Downloads.download(dictionary_url, "dictionary.txt")
end

# %%
# https://gist.github.com/jonathanBieler/de37a190d590297fc6b5d0ffee3c18dc

#=
# Port of Peter Norvig's Common Lisp program from http://norvig.com/java-lisp.html.
#
# - Julia version: 1.6.2
# - Author: Renato Athaydes
# - Date: 2021-07-24
=#
const emptyStrings = String[]

function printTranslations(num, digits, start=1, words=String[])
    if start > length(digits)
       return println(num, ": ", join(words, " "))
    end
    foundWord = false
    n = BigInt(1)
    for i in start:length(digits)
        n = n * 10 + nthDigit(digits, i)
        for word in get(dict, n, emptyStrings)
            foundWord = true
            printTranslations(num, digits, i + 1, [words; word])
        end
    end
    if !foundWord &&
        !(!isempty(words) && length(words[end]) == 1 && isdigit(words[end][begin]))
        printTranslations(num, digits, start + 1, [words; string(nthDigit(digits, start))])
    end
end

function loadDictionary(file)::Dict{BigInt, Vector{String}}
    local dict = Dict{BigInt, Vector{String}}()
    for word in eachline(file)
        push!(get!(dict, wordToNumber(word)) do; String[] end, word)
    end
    dict
end

function nthDigit(digits::String, i::Int64)::UInt
    UInt(digits[i]) - UInt('0')
end

function charToDigit(ch::Char)::UInt
    ch = lowercase(ch)
    ch == 'e' && return 0
    ch in ['j', 'n', 'q'] && return 1
    ch in ['r', 'w', 'x'] && return 2
    ch in ['d', 's', 'y'] && return 3
    ch in ['f', 't'] && return 4
    ch in ['a', 'm'] && return 5
    ch in ['c', 'i', 'v'] && return 6
    ch in ['b', 'k', 'u'] && return 7
    ch in ['l', 'o', 'p'] && return 8
    ch in ['g', 'h', 'z'] && return 9
    throw(DomainError(ch, "Not a letter"))
end

function wordToNumber(word::String)::BigInt
    n = BigInt(1)
    for ch in word
        if isletter(ch) && isascii(ch)
            n = n * 10 + charToDigit(ch)
        end
    end
    n
end


# dict = open(isempty(ARGS) ? "tests/words.txt" : ARGS[begin]) do file
#     loadDictionary(file)
# end

# open(length(ARGS) < 2 ? "tests/numbers.txt" : ARGS[begin+1]) do file
#     for num in eachline(file)
#         printTranslations(num, filter(isdigit, num))
#     end
# end

# %%
# original printTranslations added the arguments io and dict
function printTranslations_original(io, dict, num, digits, start=1, words=String[])
    if start > length(digits)
       return println(io, num, ": ", join(words, " "))
    end
    foundWord = false
    n = BigInt(1)
    for i in start:length(digits)
        n = n * 10 + nthDigit(digits, i)
        for word in get(dict, n, emptyStrings)
            foundWord = true
            printTranslations_original(io, dict, num, digits, i + 1, [words; word])
        end
    end
    if !foundWord &&
        !(!isempty(words) && length(words[end]) == 1 && isdigit(words[end][begin]))
        printTranslations_original(io, dict, num, digits, start + 1, [words; string(nthDigit(digits, start))])
    end
end

# translate numbers by printTranslations with dict and print the result to io
function translate(io::IO, printTranslations, dict, numbers)
    (num -> printTranslations(io, dict, num, filter(isdigit, num))).(numbers)
end

# %%
# generate the test data `numbers` of length n = 10^6

using Random

const _chars = ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '/', '-')
const _maxlen = 50
const _allowempty = false

function randphone(; chars=_chars, maxlen=_maxlen, allowempty=_allowempty)
    while true
        phone = randstring(chars, rand(1:maxlen))
        (allowempty || occursin(r"[0-9]", phone)) && return phone
    end
end

function randphone(n; chars=_chars, maxlen=_maxlen, allowempty=_allowempty)
    [randphone(; chars, maxlen, allowempty) for _ in 1:n]
end

n = 10^6
println("numbers = randphone($n):")
@time numbers = randphone(n)
println()

# %%
# load dictionary
println("open(loadDictionary, \"dictionary.txt\"):")
@time dict = open(loadDictionary, "dictionary.txt")
println()

# %%
# translate numbers with dict and print the result to io
io = IOBuffer()
println("translate by printTranslations_original:")
@time translate(io, printTranslations_original, dict, numbers)
println()
result_original = String(take!(io))
close(io)

# %%
# https://gist.github.com/jonathanBieler/de37a190d590297fc6b5d0ffee3c18dc

#=
# Port of Peter Norvig's Common Lisp program from http://norvig.com/java-lisp.html.
#
# - Julia version: 1.6.2
# - Author: Renato Athaydes
# - Date: 2021-07-24
=#
# const emptyStrings = String[]

function printTranslations(io, dict, num, digits, start=1, words=String[])
    if start > length(digits)
       return println(io, num, ": ", join(words, " "))
    end
    foundWord = false
    n = BigInt(1)
    for i in start:length(digits)
        Base.GMP.MPZ.mul_si!(n, 10)
        Base.GMP.MPZ.add_ui!(n, nthDigit(digits, i))
        for word in get(dict, n, emptyStrings)
            foundWord = true
            printTranslations(io, dict, num, digits, i + 1, [words; word])
        end
    end
    if !foundWord &&
        !(!isempty(words) && length(words[end]) == 1 && isdigit(words[end][begin]))
        printTranslations(io, dict, num, digits, start + 1, [words; string(nthDigit(digits, start))])
    end
end

function loadDictionary(file)::Dict{BigInt, Vector{String}}
    local dict = Dict{BigInt, Vector{String}}()
    for word in eachline(file)
        push!(get!(dict, wordToNumber(word)) do; String[] end, word)
    end
    dict
end

function nthDigit(digits::String, i::Int64)
    UInt(digits[i]) - UInt('0')
end

function charToDigit(ch::Char)
    ch = lowercase(ch)
    ch == 'e' && return 0
    ch in ('j', 'n', 'q') && return 1
    ch in ('r', 'w', 'x') && return 2
    ch in ('d', 's', 'y') && return 3
    ch in ('f', 't') && return 4
    ch in ('a', 'm') && return 5
    ch in ('c', 'i', 'v') && return 6
    ch in ('b', 'k', 'u') && return 7
    ch in ('l', 'o', 'p') && return 8
    ch in ('g', 'h', 'z') && return 9
    throw(DomainError(ch, "Not a letter"))
end

function wordToNumber(word::String)
    n = BigInt(1)
    for ch in word
        if isletter(ch) && isascii(ch)
            Base.GMP.MPZ.mul_si!(n, 10)
            Base.GMP.MPZ.add_ui!(n, charToDigit(ch))
        end
    end
    n
end

# patch in method to add integer to BigInt in-place
@eval Base.GMP.MPZ begin
    add_ui!(x::BigInt, a::BigInt, b) = (ccall((:__gmpz_add_ui, :libgmp), Cvoid, (mpz_t, mpz_t, Clong), x, a, b); x)
    add_ui!(x::BigInt, b) = add_ui!(x, x, b)
end

function main()

    path = "/Users/jbieler/Downloads/tmp/prechelt-phone-number-encoding/"

    dict = open("$(path)dictionary.txt") do file
        loadDictionary(file)
    end
    io = IOBuffer()
    open("$(path)input.txt") do file
        for num in eachline(file)
            printTranslations(io, dict, num, filter(isdigit, num))
        end
    end

end

# main() #make sure it's compiled before timing it
# @time main()

# %%
io = IOBuffer()
println("translate by printTranslations of jonathanBieler:")
@time translate(io, printTranslations, dict, numbers)
println()
result_jonathanBieler = String(take!(io))
close(io)

# %%
@show result_original == result_jonathanBieler;

# %%
