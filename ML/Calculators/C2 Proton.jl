using Plots
using Statistics
using LsqFit
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

function MLProton(plotrange)
    dir0="C:\\Users\\Drew\\Desktop\\FakeData"
    plotdir="C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\FinalPlots"

    gaugeconfigs = ["$(748 + 16*i)" for i in range(0,Integer((1420-748)/16),step=1)]
    filter!(e->e∉["956","1004","1036","1052"], gaugeconfigs)

    cd(dir0)
    filelist = readdir()
    fpconfig = Integer(length(filelist)/39)
    global datamatrix = zeros((39,64))
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
        lines=split(reldata,"\n")

        global linematrices=[] # matrix full of pieces of each line in the relevant data file
        global C2 = [] # one function C2(t) at gauge config μₐ

        rownum = 0
        for k in gaugeconfigs
            rownum += 1
            if (occursin(k,i))
                break
            end
        end

        # pushing split lines into a matrix
        for j in range(2,length(lines)-1,step=1)
            push!(linematrices,split(lines[j]))
            push!(C2, (value(linematrices[j-1][2])) ) # U-D contribution
        end
        datamatrix[rownum,:] += C2/fpconfig # row in datamatrix indicates gauge config
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
