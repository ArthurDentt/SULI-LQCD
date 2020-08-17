include("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Calculators\\C3 Axial Charge.jl")
include("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Calculators\\C3 Vector Charge.jl")
include("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Calculators\\C3 Scalar Charge.jl")
include("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Calculators\\C2 Proton.jl")

plotrange = 2:9
plotrangeproton = 1:64

println("$(MLProton(plotrangeproton))")
(VeF, VeR) = MLVector(plotrange)
println("Done with ML Vector!")
(AxF, AxR) = MLAxial(plotrange)
println("Done with ML Axial!")
(ScF, ScR) = MLScalar(plotrange)
println("Done with ML Scalar!")
(PAxF,PAxSEF) = physax(AxF,VeF)
(PAxR,PAxSER) = physax(AxR,VeR)
(PScF,PScSEF) = physax(ScF,VeF)
(PScR,PScSER) = physax(ScR,VeR)
PAxF = [mean(PAxF[:,i]) for i in range(1,length(PAxF[1,:]),step=1)]
PAxR = [mean(PAxR[:,i]) for i in range(1,length(PAxR[1,:]),step=1)]
PScF = [mean(PScF[:,i]) for i in range(1,length(PScF[1,:]),step=1)]
PScR = [mean(PScR[:,i]) for i in range(1,length(PScR[1,:]),step=1)]
plotpax(PAxR,PAxSER,PAxF,PAxSEF,plotrange)
plotpsc(PScR,PScSER,PScF,PScSEF,plotrange)
println("Done!")
