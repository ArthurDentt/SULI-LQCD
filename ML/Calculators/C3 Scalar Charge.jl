using Plots
using Statistics
using LsqFit
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

dir0="C:\\Users\\Drew\\Desktop\\FakeData"
cd(dir0)
plotdir="C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\FinalPlots"

if (isfile("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Data\\C2ProtonGGData.txt"))
    rm("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Data\\C2ProtonGGData.txt")
end

filelist = readdir()

datamatrix = zeros((length(filelist),64))

cd(dir0)
fileindex = 0
@progress (name = "$fileindex / $(length(filelist))") for i in filelist
    global fileindex += 1
    try
        global test_file=open(i,"r+");
    catch (SystemError)
        println("No file named \"$i\" exists.")
        # Un-comment to see what files are missing, I don't really care though
        continue
    end

    readuntil(test_file, "G5G3", keep = false)
    reldata = readuntil(test_file, "END_NUC3PT", keep=false)
    lines=split(reldata,"\n")

    linematrices=[] # matrix full of pieces of each line in the relevant data file
    C2 = []

    # pushing split lines into a matrix
    for i in range(6,length(lines)-1,step=1)
        push!(linematrices,split(lines[i]))
        push!(C2, (value(linematrices[i-5][2])-value(linematrices[i-5][4])) ) # U-D contribution
    end

    datamatrix[fileindex,:] = C2
    close(test_file)
end

Jackreplicates = [Jackrep(datamatrix[:,i]) for i in range(1,length(datamatrix[1,:]),step=1)]
Jackestimates = [mean(i) for i in Jackreplicates]
Stderrors = [JackSE(i) for i in Jackreplicates]

cd(plotdir)

plotrange = 1:64

scatter(plotrange[1]-1:plotrange[end]-1,Jackestimates[plotrange],marker=(:x),markercolor=(:red),
    linecolor=(:red),markerstrokecolor=(:red),yerror=Stderrors[plotrange],
    dpi=600,grid=false,frame=(:box), foreground_color_legend = nothing, background_color_legend=nothing,
    label = "170 MeV AMA", legend = ((.85,.95)), legendfontsize = 7)
xlabel!("τ");ylabel!("gₛ");title!("Scalar Charge")
savefig("Scalar Charge C3 Plot.png")
