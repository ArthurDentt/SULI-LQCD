# Plots C2(t), Finds effective mass
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

global datafile=open("C2ProtonGGData.txt","r");
datamatrix=readlines(datafile)
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

C2 = n_datamatrix[1,:]
C2SE = n_datamatrix[2,:]
Effmass = n_datamatrix[3,:]
EffmassSE = n_datamatrix[4,:]

EffectiveMass,EffectiveMassSE = massconvert(EffectiveMass,EffectiveMassSE)

cd("C:\\Users\\Drew\\github\\SULI-LQCD\\FinalPlots")

plot(1:length(C2),C2,markerstrokecolor=(:black),marker=(:x),yerror=C2SE,legend=false,dpi=600,grid=false)
xlabel!("t");ylabel!("Re(<C₂>)");title!("ProtonGG Re(<C₂>)(t)")
savefig("ProtonGG C2 Plot.png")

plot(1:length(Effmass),Effmass,markerstrokecolor=(:black),marker=(:x),legend=false,dpi=600,yerror=EffmassSE,grid=false)
xlabel!("t");ylabel!("m*");title!("ProtonGG m*(t)")
savefig("ProtonGG Emass Plot.png")

println("------------------------------------------------------")
println("Effective mass: $EffectiveMass ⁺/₋ $EffectiveMassSE MeV")

close(datafile)
println("Done")