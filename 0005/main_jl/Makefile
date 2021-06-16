#JULIA = D:/Julia-1.8.0-DEV/bin/julia
JULIA = julia

SOLS = sol_m.txt sol_t.txt sol_U.txt sol_x.txt
IMAGES = heateq.gif heateq0.png
PKG = HeatEqCalc

all:
	$(JULIA) main.jl PRECOMPILE
	$(JULIA) main.jl
	$(JULIA) animate.jl

bench:
	time julia main.jl PRECOMPILE
	time julia main.jl

clean: 
	rm -f $(SOLS) $(IMAGES)
	rm -rf $(PKG)