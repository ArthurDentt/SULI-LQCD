# Plots C3(t), Finds effective mass
using Plots
using Statistics

include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

Scale = 1
xtickvals = [i for i in range(0,9,step=1)]
ytickvals = [(1+.05*i) for i in range(0,8,step=1)]
ytick0 = ["" for i in ytickvals]
xtick0 = ["" for i in xtickvals]

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

Charges = [mean(Ratio[:,i]) for i in range(1,length(Ratio[1,:]),step=1)] # Jack replicates of charge
Charge = mean(Charges) #Jack estimator of Charge
ChargeSE = JackSE(Charges)
chisq = 0
for i in range(1,length(RatioEstimates),step=1)
    global chisq += ((Charge-RatioEstimates[i]))^2  / (RatioSE[i]) / (length(Ratio[:,1])-1)
end
println("χ²ᵥ = $chisq")

cd("C:\\Users\\Drew\\github\\SULI-LQCD\\FinalPlots")

# Defining the properly round charge and error for plot annotation
sigCharge = round(Charge,digits=3)
sigChargeSE = Integer(round(ChargeSE, digits=2)*1000)

scatter(1:length(RatioEstimates),RatioEstimates,marker=(:x),markercolor=(:red),
    linecolor=(:red),markerstrokecolor=(:red),yerror=RatioSE,
    framealpha=0,dpi=600,grid=false,frame=(:box),ylims=(ytickvals[1],ytickvals[end]),
    xlims=(xtickvals[1],xtickvals[end]),xticks = (xtickvals),
    yticks=(ytickvals), foreground_color_legend = nothing, background_color_legend=nothing, label = "170 MeV AMA",
    legend = ((.85,.95)), legendfontsize = 7)
xlabel!("τ");ylabel!("gₐ/gᵥ");title!("Physical Axial Charge")
plot!(twinx(), xmirror=:true,grid=:false,ylims=(ytickvals[1],ytickvals[end]),
    xlims=(xtickvals[1],xtickvals[end]),xticks = (xtickvals,xtick0),
    yticks=(ytickvals,ytick0))
hline!(([1.2755]),linecolor=(:blue),label="")
annotate!(0.05, 1.285, text("Experiment: 1.2732(23)",6, :left))
hline!(([1.2709]),linecolor=(:blue),label="")
hline!(([Charge + ChargeSE]),ls = :dash, lc = :black, label="")
annotate!(6.6, Charge - ChargeSE - .015, text(
    "Physical Axial Charge: $sigCharge($sigChargeSE) \n χ²ᵥ = $(round(chisq, digits=5))",6, :left))
hline!(([Charge - ChargeSE]),ls = :dash, lc = :black, label="")

savefig("Physical Axial Charge Plot.png")
println("Done")
