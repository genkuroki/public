f(x) = begin
        x2 = x ^ 2
        x * Main.muladd(x2, Main.muladd(x2, 0.008333333333333333, -0.16666666666666666), 1.0)
    end