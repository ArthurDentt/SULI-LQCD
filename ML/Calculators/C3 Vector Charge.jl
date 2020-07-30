using Plots
using Statistics
using LsqFit
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

function MLVector(plotrange)
    dirF="C:\\Users\\Drew\\Desktop\\FakeData"
    dirR="C:\\Users\\Drew\\Desktop\\RealData"
    plotdir="C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\FinalPlots"
    datadir="C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\Data"

    cd(datadir)
    global datafile=open("C2ProtonData.txt","r");
    sourcedata=readlines(datafile)
    close(datafile)
    sinkdata = split(split(split(sourcedata[1],"[")[2],"]")[1], ",")
    sinkdata = [parse(Float64,i) for i in sinkdata]

    gaugeconfigs = ["$(748 + 16*i)" for i in range(0,Integer((1420-748)/16),step=1)]
    filter!(e->e∉["956","1004","1036","1052"], gaugeconfigs)

    # FAKE DATA PROCESS
    cd(dirF)
    filelist = readdir()
    fpconfig = Integer(length(filelist)/length(gaugeconfigs))
    global datamatrixF = zeros((length(gaugeconfigs),64))
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

        readuntil(test_file, "G4", keep = false)
        reldata = readuntil(test_file, "END_NUC3PT", keep=false)
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
        for j in range(6,length(lines)-1,step=1)
            push!(linematrices,split(lines[j]))
            push!(C2, (value(linematrices[j-5][2])-value(linematrices[j-5][4])) ) # U-D contribution
        end
        datamatrixF[rownum,:] += C2/(fpconfig*abs(sinkdata[rownum])) # row in datamatrix indicates gauge config
        close(test_file)
    end

    # REAL DATA PROCESS
    cd(dirR)
    filelist = readdir()
    fpconfig = Integer(length(filelist)/length(gaugeconfigs))
    global datamatrixR = zeros((length(gaugeconfigs),64))
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

        readuntil(test_file, "G4", keep = false)
        reldata = readuntil(test_file, "END_NUC3PT", keep=false)
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
        for j in range(6,length(lines)-1,step=1)
            push!(linematrices,split(lines[j]))
            push!(C2, (value(linematrices[j-5][2])-value(linematrices[j-5][4])) ) # U-D contribution
        end
        datamatrixR[rownum,:] += C2/(fpconfig*abs(sinkdata[rownum])) # row in datamatrix indicates gauge config
        close(test_file)
    end

    JackreplicatesF = [Jackrep(datamatrixF[:,i]) for i in range(1,length(datamatrixF[1,:]),step=1)]
    global JackestimatesF = [mean(i) for i in JackreplicatesF]
    StderrorsF = [JackSE(i) for i in JackreplicatesF]

    JackreplicatesR = [Jackrep(datamatrixR[:,i]) for i in range(1,length(datamatrixR[1,:]),step=1)]
    JackestimatesR = [mean(i) for i in JackreplicatesR]
    StderrorsR = [JackSE(i) for i in JackreplicatesR]

    cd(plotdir)

    scatter(plotrange[1]-1:plotrange[end]-1,JackestimatesR[plotrange],marker=(:x),markercolor=(:purple),
        linecolor=(:purple),markerstrokecolor=(:purple),yerror=StderrorsR[plotrange],
        dpi=600,grid=false,frame=(:box), foreground_color_legend = nothing, background_color_legend=nothing,
        label = "REAL (170 MeV AMA)", legend = ((.75,.95)), legendfontsize = 7)
    scatter!(plotrange[1]-1:plotrange[end]-1,JackestimatesF[plotrange],marker=(:x),markercolor=(:green),
        linecolor=(:green),markerstrokecolor=(:green),yerror=StderrorsF[plotrange],
        label = "ML (170 MeV AMA)")
    xlabel!("τ");ylabel!("gᵥ");title!("Vector Charge")
    savefig("Vector Charge C3 Plot.png")

    return("Done with ML Vector!")
end
