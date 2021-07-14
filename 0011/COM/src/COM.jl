"""
`COM.prettyprint(x)` prints nicely the object `x` for which `COM.prettify(x)` method is defined.
"""
module COM
prettify(x) = sprint(io -> show(io, "text/plain", x))
prettyprint(io::IO, x) = print(io, prettify(x))
prettyprint(x) = prettyprint(stdout, x)
end
