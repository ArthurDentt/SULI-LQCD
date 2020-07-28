using Plots
using Statistics
using LsqFit
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

function MLProton(plotrange)
    dir0="C:\\Users\\Drew\\Desktop\\RealData"
    cd(dir0)
    plotdir="C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\FinalPlots"

    filelist = readdir()

    datamatrix = zeros((length(filelist),64))

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

        readuntil(test_file, "NUC_G5C_PP", keep = false)
        reldata = readuntil(test_file, "ENDPROP", keep=false)
        global lines=split(reldata,"\n")

        global linematrices=[] # matrix full of pieces of each line in the relevant data file
        global C2 = []

        # pushing split lines into a matrix
        for i in range(2,length(lines)-1,step=1)
            push!(linematrices,split(lines[i]))
            push!(C2, value(linematrices[i-1][2])) # positive parity real part only
        end

        datamatrix[fileindex,:] = C2
        close(test_file)
    end

    Jackreplicates = [Jackrep(datamatrix[:,i]) for i in range(1,length(datamatrix[1,:]),step=1)]
    Jackestimates = [mean(i) for i in Jackreplicates]
    Stderrors = [JackSE(i) for i in Jackreplicates]

    cd(plotdir)

    plot(plotrange[1]-1:plotrange[end]-1,Jackestimates[plotrange],marker=(:x),
        markerstrokecolor=(:black),yerror=Stderrors[plotrange],
        dpi=600,grid=false,frame=(:box), foreground_color_legend = nothing, background_color_legend=nothing,
        label = "170 MeV AMA", legend = ((.85,.95)), legendfontsize = 7)
    xlabel!("t");ylabel!("Re(C2)");title!("C2 Proton")
    savefig("C2 Proton Plot.png")

    return("Done with ML Proton!")
end
