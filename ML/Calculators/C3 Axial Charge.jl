using Plots
using Statistics
using LsqFit
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

function MLAxial(plotrange)
    dir0="C:\\Users\\Drew\\Desktop\\RealData"
    cd(dir0)
    plotdir="C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\FinalPlots"

    if (isfile("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Data\\C2ProtonGGData.txt"))
        rm("C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Data\\C2ProtonGGData.txt")
    end

    filelist = readdir()

    global datamatrix = zeros((length(filelist),64))

    cd(dir0)
    fileindex = 0
    @progress (name = "$fileindex / $(length(filelist))") for i in filelist
        fileindex += 1
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

        global linematrices=[] # matrix full of pieces of each line in the relevant data file
        global C2 = [] # one function C2(t) at gauge config μₐ

        # pushing split lines into a matrix
        for i in range(6,length(lines)-1,step=1)
            push!(linematrices,split(lines[i]))
            push!(C2, (value(linematrices[i-5][2])-value(linematrices[i-5][4])) ) # U-D contribution
        end

        datamatrix[fileindex,:] = C2 # row in datamatrix indicates gauge config
        close(test_file)
    end

    Jackreplicates = [Jackrep(datamatrix[:,i]) for i in range(1,length(datamatrix[1,:]),step=1)]
    Jackestimates = [mean(i) for i in Jackreplicates]
    Stderrors = [JackSE(i) for i in Jackreplicates]

    cd(plotdir)

    scatter(plotrange[1]-1:plotrange[end]-1,Jackestimates[plotrange],marker=(:x),markercolor=(:red),
        linecolor=(:red),markerstrokecolor=(:red),yerror=Stderrors[plotrange],
        dpi=600,grid=false,frame=(:box), foreground_color_legend = nothing, background_color_legend=nothing,
        label = "170 MeV AMA", legend = ((.85,.95)), legendfontsize = 7)
    xlabel!("τ");ylabel!("gₐ");title!("Axial Charge")
    savefig("Axial Charge C3 Plot.png")

    return("Done with ML Axial!")
end
