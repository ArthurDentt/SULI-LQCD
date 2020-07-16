# Plots C3(t), Finds effective mass
using Plots
using Statistics

Scale = sqrt(2)
plotrange = 1:9

xtickvals = [i for i in range(0,9,step=1)]
ytickvals = [(0 +.5*i) for i in range(0,4,step=1)]
ytick0 = ["" for i in range(1,length(ytickvals),step=1)]
xtick0 = ["" for i in range(1,length(xtickvals),step=1)]

xtickvals2 = [i for i in range(0,9,step=1)]
ytickvals2 = [(0.6 +.2*i) for i in range(0,6,step=1)]
ytick02 = ["" for i in range(1,length(ytickvals2),step=1)]
xtick02 = ["" for i in range(1,length(xtickvals),step=1)]

println("Starting output...")

cd("C:\\Users\\Drew\\github\\SULI-LQCD\\Data")

global datafile=open("C3ScalarChargeData.txt","r");
datamatrix=readlines(datafile)
close(datafile)
datamatrix = [split(split(split(datamatrix[i],"[")[2],"]")[1], ",") for i in range(1,4,step=1)]

# Take incoming data and reformat it into a matrix where each row is a datavector
# rows: Data 1, Error 1, Data 2, Error 2, ... etc.
n_datamatrix  = zeros((length(datamatrix),length(datamatrix[1])))
for i in range(1,length(datamatrix),step=1)
    datavector = zeros(length(datamatrix[i]))
    for j in range(1,length(datamatrix[i]),step=1)
        datavector[j] = parse(Float64,datamatrix[i][j])
    end
    n_datamatrix[i,1:length(datavector)]=datavector
end
C3 = []
C3SE = []

for i in range(1,length(n_datamatrix[1,:]),step=1)
    if n_datamatrix[1,i] < 50
        push!(C3,n_datamatrix[1,i])
        push!(C3SE,n_datamatrix[2,i])
    end
end
physgs = n_datamatrix[3,1:8]
physgsSE = n_datamatrix[4,1:8]
cd("C:\\Users\\Drew\\github\\SULI-LQCD\\FinalPlots")

scatter(0:length(plotrange)-1,C3[plotrange]*Scale,marker=(:x),markercolor=(:red),
    linecolor=(:red),markerstrokecolor=(:red),yerror=C3SE[plotrange]*Scale,
    dpi=600,grid=false,xlims=(xtickvals[1],xtickvals[end]),
    ylims=(ytickvals[1],ytickvals[end]),xticks=xtickvals,yticks=ytickvals,frame=(:box),
    foreground_color_legend = nothing, label = "170 MeV AMA",
    legend = ((.85,.95)), legendfontsize = 7)
xlabel!("τ");ylabel!("gₛ");title!("Scalar Charge")
plot!(twinx(), xmirror=:true,grid=:false,ylims=(ytickvals[1],ytickvals[end]),
    xlims=(xtickvals[1],xtickvals[end]),xticks = (xtickvals,xtick0),
    yticks=(ytickvals,ytick0))
savefig("Scalar Charge C3 Plot.png")

scatter(1:8,physgs*Scale,marker=(:x),markercolor=(:red),
    linecolor=(:red),markerstrokecolor=(:red),yerror=physgsSE*Scale,
    dpi=600,grid=false,xlims=(xtickvals2[1],xtickvals2[end]),
    ylims=(ytickvals2[1],ytickvals2[end]),xticks=xtickvals2,yticks=ytickvals2,frame=(:box),
    foreground_color_legend = nothing, label = "170 MeV AMA",
    legend = ((.85,.95)), legendfontsize = 7)
xlabel!("τ");ylabel!("gₛZₛ");title!("Physical Scalar Charge")
plot!(twinx(), xmirror=:true,grid=:false,ylims=(ytickvals2[1],ytickvals2[end]),
    xlims=(xtickvals2[1],xtickvals2[end]),xticks = (xtickvals2,xtick02),
    yticks=(ytickvals2,ytick02))
savefig("Physical Scalar Charge Plot.png")

println("------------------------------------------------------")
close(datafile)
println("Done")
