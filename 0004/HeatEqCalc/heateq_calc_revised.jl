using Pkg; Pkg.activate(".")
using HeatEqCalc
sol = HeatEqCalc.calc_sol()
HeatEqCalc.save_sol(sol)
