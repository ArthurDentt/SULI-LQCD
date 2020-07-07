# Plots C3(t), Finds effective mass
using Plots
using Statistics

# Converts from unitless to MeV
function massconvert(mass,massSE)
    invspacing = 1.378
    invspacingSE = .007
    n_mass = mass*invspacing
    n_massSE = sqrt( (mass*invspacingSE)^2 + (invspacing*massSE)^2 )
    return([1000*n_mass,1000*n_massSE])
end

println("Starting output...")

cd("C:\\Users\\Drew\\github\\SULI-LQCD")

C2T=188167.6457778
C2TSE=982.7230311587001

global datafile=open("C3VectorChargeData.txt","r");
datamatrix=readlines(datafile)
close(datafile)
EffectiveMass = parse(Float64,datamatrix[5])
EffectiveMassSE = parse(Float64,datamatrix[6])
datamatrix = [split(split(split(datamatrix[i],"[")[2],"]")[1], ",") for i in range(1,4,step=1)]

# Take incoming data and reformat it into a matrix where each row is a datavector
# rows: Data 1, Error 1, Data 2, Error 2, ... etc.
n_datamatrix  = zeros((length(datamatrix),length(datamatrix[1])))
for i in range(1,length(datamatrix),step=1)
    datavector = zeros(length(datamatrix[i]))
    for j in range(1,length(datamatrix[i]),step=1)
        datavector[j] = parse(Float64,datamatrix[i][j])
    end
    n_datamatrix[i,:]=datavector
end
C3 = []
C3SE = []
n_datamatrix[1,:]
for i in range(1,length(n_datamatrix[1,:]),step=1)
    if n_datamatrix[1,i] < 10^7
        push!(C3,n_datamatrix[1,i])
        push!(C3SE,n_datamatrix[2,i])
    end
end

Effmass = n_datamatrix[3,:]
EffmassSE = n_datamatrix[4,:]

EffectiveMass,EffectiveMassSE = massconvert(EffectiveMass,EffectiveMassSE)
cd("C:\\Users\\Drew\\github\\SULI-LQCD\\FinalPlots")

plot(1:length(C3),C3,markerstrokecolor=(:black),marker=(:circle),yerror=C3SE,legend=false,dpi=600,grid=false)
xlabel!("τ (τ=10 has been removed)");ylabel!("C₃/C₂(T)");title!("Vector Charge C₃(τ,T)/C₂(T)")
savefig("Vector Charge C3 Plot.png")

plot(1:length(Effmass),Effmass,markerstrokecolor=(:black),marker=(:circle),legend=false,dpi=600,yerror=EffmassSE,grid=false)
xlabel!("τ");ylabel!("m*");title!("Vector Charge m*(τ,T)")
savefig("Vector Charge Emass Plot.png")

println("------------------------------------------------------")
println("Effective mass: $EffectiveMass ⁺/₋ $EffectiveMassSE MeV (This has no meaning right now)")

close(datafile)
println("Done")
