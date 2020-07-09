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
Scale = sqrt(2)
plotrange = 1:9

xtickvals = [i for i in range(0,9,step=1)]
ytickvals = [(0 +.5*i) for i in range(0,4,step=1)]
ytick0 = ["" for i in range(1,length(ytickvals),step=1)]
ytick0 = ["" for i in range(1,length(xtickvals),step=1)]

println("Starting output...")

cd("C:\\Users\\Drew\\github\\SULI-LQCD")

global datafile=open("C3ScalarChargeData.txt","r");
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
    if n_datamatrix[1,i] < 50
        push!(C3,n_datamatrix[1,i])
        push!(C3SE,n_datamatrix[2,i])
    end
end

Effmass = n_datamatrix[3,:]
EffmassSE = n_datamatrix[4,:]

EffectiveMass,EffectiveMassSE = massconvert(EffectiveMass,EffectiveMassSE)
cd("C:\\Users\\Drew\\github\\SULI-LQCD\\FinalPlots")

scatter(0:length(plotrange)-1,C3[plotrange]*Scale,marker=(:x),markercolor=(:red),
    linecolor=(:red),markerstrokecolor=(:red),yerror=C3SE[plotrange]*Scale,
    legend=false,dpi=600,grid=false,xlims=(xtickvals[1],xtickvals[end]),
    ylims=(ytickvals[1],ytickvals[end]),xticks=xtickvals,yticks=ytickvals,frame=(:box))
xlabel!("τ");ylabel!("gₛ");title!("Scalar Charge")
plot!(twinx(), xmirror=:true,grid=:false,ylims=(ytickvals[1],ytickvals[end]),
    xlims=(xtickvals[1],xtickvals[end]),xticks = (xtickvals,xtick0),
    yticks=(ytickvals,ytick0))
savefig("Scalar Charge C3 Plot.png")

println("------------------------------------------------------")
println("Effective mass: $EffectiveMass ⁺/₋ $EffectiveMassSE MeV (This has no meaning right now)")

close(datafile)
println("Done")
