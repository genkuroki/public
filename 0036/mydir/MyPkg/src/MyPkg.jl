module MyPkg

using Base.Threads

greet() = print("Hello My Package!")

function mcpi(N)
    tmp = zeros(Int, nthreads())
    @threads for _ in 1:N
        tmp[threadid()] += rand()^2 + rand()^2 ≤ 1
    end
    4sum(tmp)/N
end

end # module MyPkg
