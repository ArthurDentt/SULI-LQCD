# Plots C2(t), Finds effective mass
using Plots
using Statistics

println("Starting output...")

cd("C:\\Users\\Drew\\github\\SULI-LQCD")

global datafile=open("ProtonC2Data.txt","r");
datamatrix=readlines(datafile)
datamatrix = [split(split(split(split(datamatrix[i],"Any")[2],"[")[2],"]")[1], ",") for i in range(1,length(datamatrix),step=1)]

slices = [4:9,11:16,21:26,27:32,35:40,43:48,51:56] # Slices for probing Emass

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

global Tindex = 0
cd("C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex")
rdir = [i for i in readdir()]
global k=1
@progress (name="Plotting...") for i in range(1,length(n_datamatrix[:,1]),step=2)
    cd("C:\\Users\\Drew\\Desktop\\C2Graphs\\Proton_PP\\RealEC2(t)")
    reC3plot=plot(0:length(n_datamatrix[i,:])-1,n_datamatrix[i,:],marker=(:circle),legend=false,yerror=n_datamatrix[i+1,:])
    xlabel!("t");ylabel!("Re(<C₂>)");title!("Proton_PP Re(<C₂>)(t)")
    # Displaying plots takes a while, so I don't do it
    #display(meC3plot)
    suffix = rdir[k]
    savefig("T=$(Tindex), $(suffix).png")
    #println("plotted C:\\Users\\Drew\\Desktop\\C3Graphs\\Scalar Charge\\RealEC3(t)\\T=$(Tindex), $(suffix).png")
    # Not needed due to loading bar :)
    if (k==16)
        global k=1
        global Tindex = Tindex + 8
    else
        global k=k+1
    end
end
close(datafile)
println("")
println("Done")
println("** Plotted with NO Error bars **")
println("")
