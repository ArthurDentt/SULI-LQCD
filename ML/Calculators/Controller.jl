include("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Calculators\\C3 Axial Charge.jl")
include("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Calculators\\C3 Vector Charge.jl")
include("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Calculators\\C3 Scalar Charge.jl")
include("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Calculators\\C2 Proton.jl")

plotrange = 1:9
plotrangeproton = 1:64

println("$(MLProton(plotrangeproton))")
println("$(MLAxial(plotrange))")
println("$(MLScalar(plotrange))")
println("$(MLVector(plotrange))")
println("Done!")
