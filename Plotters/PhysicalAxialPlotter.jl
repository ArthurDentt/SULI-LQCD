# Plots C3(t), Finds effective mass
using Plots
using Statistics
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

Scale = 1

println("Starting output...")

cd("C:\\Users\\Drew\\github\\SULI-LQCD\\Data")

#get gₐ,gᵥ data from files
global datafile=open("gARatioData.txt","r");
gA=readlines(datafile)
close(datafile)
gA = [split(split(split(gA[i],"[")[2],"]")[1], ",") for i in range(1,8,step=1)]
global datafile=open("gVRatioData.txt","r");
gV=readlines(datafile)
close(datafile)
gV = [split(split(split(gV[i],"[")[2],"]")[1], ",") for i in range(1,8,step=1)]

# reformatting replicates into ratio matrix -> still need to do JackSE on this
Ratio = zeros((length(gA),length(gA[1])))
for i in range(1,length(Ratio[:,1]),step=1)
    for j in range(1,length(Ratio[1,:]),step=1)
        Ratio[i,j] = parse(Float64,gA[i][j])/parse(Float64,gV[i][j])
    end
end
RatioEstimates = [mean(Ratio[i,:]) for i in range(1,length(Ratio[:,1]),step=1)]
RatioSE = [JackSE(Ratio[i,:]) for i in range(1,length(Ratio[:,1]),step=1)]
cd("C:\\Users\\Drew\\github\\SULI-LQCD\\FinalPlots")

scatter(1:length(RatioEstimates),RatioEstimates,marker=(:x),markercolor=(:red),
    linecolor=(:red),markerstrokecolor=(:red),yerror=RatioSE,
    legend=false,dpi=600,grid=false,frame=(:box))
xlabel!("τ");ylabel!("gₐ/gᵥ");title!("Physical Axial Charge")

savefig("Physical Axial Charge Plot.png")
println("Done")
