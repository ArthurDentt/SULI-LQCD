# Plots C3(t), Finds effective mass
using Plots
using Statistics

println("Starting output...")

cd("C:\\Users\\Drew\\github\\SULI-LQCD")

global datafile=open("C3ScalarChargeData.txt","r");
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

C3 = n_datamatrix[1,:]
C3SE = n_datamatrix[2,:]
Effmass = n_datamatrix[3,:]
EffmassSE = n_datamatrix[4,:]

cd("C:\\Users\\Drew\\github\\SULI-LQCD\\FinalPlots")

plot(1:length(C3),C3,markerstrokecolor=(:black),marker=(:circle),yerror=C3SE,legend=false,dpi=600,grid=false)
xlabel!("t");ylabel!("Re(<C₃>)");title!("Scalar Charge Re(<C₃>)(t)")
savefig("Scalar Charge C3 Plot.png")

plot(1:length(Effmass),Effmass,markerstrokecolor=(:black),marker=(:circle),legend=false,dpi=600,yerror=EffmassSE,grid=false)
xlabel!("t");ylabel!("m*");title!("Scalar Charge m*(t)")
savefig("Scalar Charge Emass Plot.png")

println("------------------------------------------------------")
println("Effective mass: $EffectiveMass ⁺/₋ $EffectiveMassSE")

close(datafile)
println("Done")
