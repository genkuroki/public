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

# %% [markdown]
# For REPL, see https://docs.julialang.org/en/v1/stdlib/REPL/#The-Julian-mode

# %%
x = 1.234567890123456789

# %%
const COMPACT = Ref(true)

# copy from
# https://github.com/JuliaLang/IJulia.jl/blob/cc2a9bf61a2515596b177339f9a3514de8c38573/src/inline.jl#L29
#
@eval IJulia InlineIOContext(io, KVs::Pair...) = IOContext(
    io,
    :limit=>true, :color=>true, :jupyter=>true,
    :compact=>Main.COMPACT[],
    KVs...
)

# %%
COMPACT[] = true
x

# %%
COMPACT[] = false
x

# %%
buf = Vector{UInt8}(undef, 20)
pos = Base.Ryu.writeshortest(buf, 1, x, false, false, true, -1, 
    UInt8('f'), false, UInt8('.'), false, false)
sprint(write, resize!(buf, pos - 1))

# %%
buf = Vector{UInt8}(undef, 20)
pos = Base.Ryu.writeshortest(buf, 1, x, false, false, true, -1, 
    UInt8('f'), false, UInt8('.'), false, true)
sprint(write, resize!(buf, pos - 1))

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/82ae530db62f364853b355bcd64b28c71445a2c4/base/ryu/Ryu.jl#L111

# %%
@which Base.Ryu.writeshortest(buf, 1, x, false, false, true, -1, 
    UInt8('f'), false, UInt8('.'), false, false)

# %%
Base.Ryu.reduce_shortest(x, 999_999)

# %%
@which Base.Ryu.reduce_shortest(x, 999_999)

# %%
const MAXSIGNIF = Ref(999_999)

# copy from 
# https://github.com/JuliaLang/julia/blob/82ae530db62f364853b355bcd64b28c71445a2c4/base/ryu/shortest.jl#L227
# 
@eval Base.Ryu function writeshortest(buf::Vector{UInt8}, pos, x::T,
                       plus=false, space=false, hash=true,
                       precision=-1, expchar=UInt8('e'), padexp=false, decchar=UInt8('.'),
                       typed=false, compact=false) where {T}
    @assert 0 < pos <= length(buf)
    # special cases
    if x == 0
        if typed && x isa Float16
            buf[pos] = UInt8('F')
            buf[pos + 1] = UInt8('l')
            buf[pos + 2] = UInt8('o')
            buf[pos + 3] = UInt8('a')
            buf[pos + 4] = UInt8('t')
            buf[pos + 5] = UInt8('1')
            buf[pos + 6] = UInt8('6')
            buf[pos + 7] = UInt8('(')
            pos += 8
        end
        pos = append_sign(x, plus, space, buf, pos)
        buf[pos] = UInt8('0')
        pos += 1
        if hash
            buf[pos] = decchar
            pos += 1
        end
        if precision == -1
            buf[pos] = UInt8('0')
            pos += 1
            if typed && x isa Float32
                buf[pos] = UInt8('f')
                buf[pos + 1] = UInt8('0')
                pos += 2
            end
            if typed && x isa Float16
                buf[pos] = UInt8(')')
                pos += 1
            end
            return pos
        end
        while hash && precision > 1
            buf[pos] = UInt8('0')
            pos += 1
            precision -= 1
        end
        if typed && x isa Float32
            buf[pos] = UInt8('f')
            buf[pos + 1] = UInt8('0')
            pos += 2
        end
        if typed && x isa Float16
            buf[pos] = UInt8(')')
            pos += 1
        end
        return pos
    elseif isnan(x)
        pos = append_sign(x, plus, space, buf, pos)
        buf[pos] = UInt8('N')
        buf[pos + 1] = UInt8('a')
        buf[pos + 2] = UInt8('N')
        if typed
            if x isa Float32
                buf[pos + 3] = UInt8('3')
                buf[pos + 4] = UInt8('2')
            elseif x isa Float16
                buf[pos + 3] = UInt8('1')
                buf[pos + 4] = UInt8('6')
            end
        end
        return pos + 3 + (typed && x isa Union{Float32, Float16} ? 2 : 0)
    elseif !isfinite(x)
        pos = append_sign(x, plus, space, buf, pos)
        buf[pos] = UInt8('I')
        buf[pos + 1] = UInt8('n')
        buf[pos + 2] = UInt8('f')
        if typed
            if x isa Float32
                buf[pos + 3] = UInt8('3')
                buf[pos + 4] = UInt8('2')
            elseif x isa Float16
                buf[pos + 3] = UInt8('1')
                buf[pos + 4] = UInt8('6')
            end
        end
        return pos + 3 + (typed && x isa Union{Float32, Float16} ? 2 : 0)
    end

    #output, nexp = reduce_shortest(x, compact ? 999_999 : nothing)
    output, nexp = reduce_shortest(x, compact ? Main.MAXSIGNIF[] : nothing)

    if typed && x isa Float16
        buf[pos] = UInt8('F')
        buf[pos + 1] = UInt8('l')
        buf[pos + 2] = UInt8('o')
        buf[pos + 3] = UInt8('a')
        buf[pos + 4] = UInt8('t')
        buf[pos + 5] = UInt8('1')
        buf[pos + 6] = UInt8('6')
        buf[pos + 7] = UInt8('(')
        pos += 8
    end
    pos = append_sign(x, plus, space, buf, pos)

    olength = decimallength(output)
    exp_form = true
    pt = nexp + olength
    if -4 < pt <= (precision == -1 ? (T == Float16 ? 3 : 6) : precision) &&
        !(pt >= olength && abs(mod(x + 0.05, 10^(pt - olength)) - 0.05) > 0.05)
        exp_form = false
        if pt <= 0
            buf[pos] = UInt8('0')
            pos += 1
            buf[pos] = decchar
            pos += 1
            for _ = 1:abs(pt)
                buf[pos] = UInt8('0')
                pos += 1
            end
            # elseif pt >= olength
            # nothing to do at this point
            # else
            # nothing to do at this point
        end
    else
        pos += 1
    end
    i = 0
    ptr = pointer(buf)
    ptr2 = pointer(DIGIT_TABLE)
    if (output >> 32) != 0
        q = output ÷ 100000000
        output2 = (output % UInt32) - UInt32(100000000) * (q % UInt32)
        output = q

        c = output2 % UInt32(10000)
        output2 = div(output2, UInt32(10000))
        d = output2 % UInt32(10000)
        c0 = (c % 100) << 1
        c1 = (c ÷ 100) << 1
        d0 = (d % 100) << 1
        d1 = (d ÷ 100) << 1
        memcpy(ptr, pos + olength - 2, ptr2, c0 + 1, 2)
        memcpy(ptr, pos + olength - 4, ptr2, c1 + 1, 2)
        memcpy(ptr, pos + olength - 6, ptr2, d0 + 1, 2)
        memcpy(ptr, pos + olength - 8, ptr2, d1 + 1, 2)
        i += 8
    end
    output2 = output % UInt32
    while output2 >= 10000
        c = output2 % UInt32(10000)
        output2 = div(output2, UInt32(10000))
        c0 = (c % 100) << 1
        c1 = (c ÷ 100) << 1
        memcpy(ptr, pos + olength - i - 2, ptr2, c0 + 1, 2)
        memcpy(ptr, pos + olength - i - 4, ptr2, c1 + 1, 2)
        i += 4
    end
    if output2 >= 100
        c = (output2 % UInt32(100)) << 1
        output2 = div(output2, UInt32(100))
        memcpy(ptr, pos + olength - i - 2, ptr2, c + 1, 2)
        i += 2
    end
    if output2 >= 10
        c = output2 << 1
        buf[pos + 1] = DIGIT_TABLE[c + 2]
        buf[pos - exp_form] = DIGIT_TABLE[c + 1]
    else
        buf[pos - exp_form] = UInt8('0') + (output2 % UInt8)
    end

    if !exp_form
        if pt <= 0
            pos += olength
            precision -= olength
            while hash && precision > 0
                buf[pos] = UInt8('0')
                pos += 1
                precision -= 1
            end
        elseif pt >= olength
            pos += olength
            precision -= olength
            for _ = 1:nexp
                buf[pos] = UInt8('0')
                pos += 1
                precision -= 1
            end
            if hash
                buf[pos] = decchar
                pos += 1
                if precision < 0
                    buf[pos] = UInt8('0')
                    pos += 1
                end
                while precision > 0
                    buf[pos] = UInt8('0')
                    pos += 1
                    precision -= 1
                end
            end
        else
            pointoff = olength - abs(nexp)
            memmove(ptr, pos + pointoff + 1, ptr, pos + pointoff, olength - pointoff + 1)
            buf[pos + pointoff] = decchar
            pos += olength + 1
            precision -= olength
            while hash && precision > 0
                buf[pos] = UInt8('0')
                pos += 1
                precision -= 1
            end
        end
        if typed && x isa Float32
            buf[pos] = UInt8('f')
            buf[pos + 1] = UInt8('0')
            pos += 2
        end
    else
        if olength > 1 || hash
            buf[pos] = decchar
            pos += olength
            precision -= olength
        end
        if hash && olength == 1
            buf[pos] = UInt8('0')
            pos += 1
        end
        while hash && precision > 0
            buf[pos] = UInt8('0')
            pos += 1
            precision -= 1
        end

        buf[pos] = expchar
        pos += 1
        exp2 = nexp + olength - 1
        if exp2 < 0
            buf[pos] = UInt8('-')
            pos += 1
            exp2 = -exp2
        elseif padexp
            buf[pos] = UInt8('+')
            pos += 1
        end

        if exp2 >= 100
            c = exp2 % 10
            memcpy(ptr, pos, ptr2, 2 * div(exp2, 10) + 1, 2)
            buf[pos + 2] = UInt8('0') + (c % UInt8)
            pos += 3
        elseif exp2 >= 10
            memcpy(ptr, pos, ptr2, 2 * exp2 + 1, 2)
            pos += 2
        else
            if padexp
                buf[pos] = UInt8('0')
                pos += 1
            end
            buf[pos] = UInt8('0') + (exp2 % UInt8)
            pos += 1
        end
    end
    if typed && x isa Float16
        buf[pos] = UInt8(')')
        pos += 1
    end

    return pos
end

# %%
A = rand(2, 2)

# %%
MAXSIGNIF[] = 99
A

# %%
MAXSIGNIF[] = 999_999
A

# %%
MAXSIGNIF[] = 999_999_999
A

# %%
COMPACT[] = true
MAXSIGNIF[] = 99
x

# %%
COMPACT[] = true
MAXSIGNIF[] = 999
x

# %%
COMPACT[] = true
MAXSIGNIF[] = 999_999
x

# %%
COMPACT[] = true
MAXSIGNIF[] = 999_999_999
x

# %%
