g(x) = begin
        begin
            while false
            end
            local var"#48#stats" = Base.gc_num()
            local var"#51#compile_elapsedtime" = Base.cumulative_compile_time_ns_before()
            local var"#50#elapsedtime" = Base.time_ns()
            local var"#49#val" = begin
                        sleep(1)
                        sin(x)
                    end
            var"#50#elapsedtime" = Base.time_ns() - var"#50#elapsedtime"
            var"#51#compile_elapsedtime" = Base.cumulative_compile_time_ns_after() - var"#51#compile_elapsedtime"
            local var"#52#diff" = Base.GC_Diff(Base.gc_num(), var"#48#stats")
            Base.time_print(var"#50#elapsedtime", (var"#52#diff").allocd, (var"#52#diff").total_time, Base.gc_alloc_count(var"#52#diff"), var"#51#compile_elapsedtime", true)
            var"#49#val"
        end
    end