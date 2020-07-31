using Plots
using Statistics
using LsqFit
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

function MLAxial(plotrange)
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

        readuntil(test_file, "G5G3", keep = false)
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
        datamatrixF[rownum,:] += C2*sqrt(2)/(fpconfig*abs(sinkdata[rownum])*3.2) # row in datamatrix indicates gauge config
        close(test_file)
    end

    # REAL DATA PROCESS
    cd(dirR)
    filelist = readdir()
    fpconfig = Integer(length(filelist)/39)
    global datamatrixR = zeros((39,64))
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
        datamatrixR[rownum,:] += C2*sqrt(2)/(fpconfig*abs(sinkdata[rownum])*3.2) # row in datamatrix indicates gauge config
        close(test_file)
    end

    JackreplicatesF = [Jackrep(datamatrixF[:,i]) for i in range(1,length(datamatrixF[1,:]),step=1)]
    global JackestimatesF = [mean(i) for i in JackreplicatesF]
    StderrorsF = [JackSE(i) for i in JackreplicatesF]

    JackreplicatesR = [Jackrep(datamatrixR[:,i]) for i in range(1,length(datamatrixR[1,:]),step=1)]
    JackestimatesR = [mean(i) for i in JackreplicatesR]
    StderrorsR = [JackSE(i) for i in JackreplicatesR]

    # Finding fake and real Charge Jackknife replicates
    plateaulength = 6
    ChargerepsF = zeros((39,plateaulength))
    ChargerepsR = zeros((39,plateaulength))
    for i in range(3,8,step=1) # 2->7
        ChargerepsF[:,i-2] = JackreplicatesF[i]
        ChargerepsR[:,i-2] = JackreplicatesR[i]
    end
    ChargerepsF = [mean(ChargerepsF[i,:]) for i in range(1,39,step=1)]
    ChargerepsR = [mean(ChargerepsR[i,:]) for i in range(1,39,step=1)]

    # Finding fake and real Charge estimates
    ChargeF = mean(ChargerepsF)
    ChargeR = mean(ChargerepsR)

    # Finding fake and real Charge SE
    ChargeSEF = JackSE(ChargerepsF)
    ChargeSER = JackSE(ChargerepsR)

    # Rounding fake and real charges and SE's to 3 digits for plot viewing
    sigChargeF = round(ChargeF,digits=3)
    sigChargeR = round(ChargeR,digits=3)
    sigChargeSEF = round(ChargeSEF,digits=3)
    sigChargeSER = round(ChargeSER,digits=3)

    # Calculating fake and real chisquared values
    chisqF = 0
    chisqR = 0
    for i in range(1,plateaulength,step=1)
        chisqF += ((ChargeF-JackestimatesF[i+2]))^2  / (StderrorsF[i+2]) / (plateaulength-1)
        chisqR += ((ChargeR-JackestimatesR[i+2]))^2  / (StderrorsR[i+2]) / (plateaulength-1)
    end
    chisqF = round(chisqF,digits=3)
    chisqR = round(chisqR,digits=3)

    cd(plotdir)

    # Setting plot boundaries
    xtickvals = [i for i in range(0,9,step=1)]
    ytickvals = [(1.45 +.05*i) for i in range(0,10,step=1)]
    ytick0 = ["" for i in range(1,length(ytickvals),step=1)]
    xtick0 = ["" for i in range(1,length(xtickvals),step=1)]

    # Plotting ML and Real comparison
    scatter(plotrange[1]-1:plotrange[end]-1,JackestimatesR[plotrange],marker=(:x),markercolor=(:purple),
        linecolor=(:purple),markerstrokecolor=(:purple),yerror=StderrorsR[plotrange],
        dpi=600,grid=false,xlims=(xtickvals[1],xtickvals[end]),
        ylims=(ytickvals[1],ytickvals[end]),xticks=xtickvals,yticks=ytickvals,frame=(:box),
        foreground_color_legend = nothing, background_color_legend=nothing, label = "REAL",
        legend = ((.85,.95)), legendfontsize = 7)
    scatter!(plotrange[1]-1:plotrange[end]-1,JackestimatesF[plotrange],marker=(:x),markercolor=(:green),
        linecolor=(:green),markerstrokecolor=(:green),yerror=StderrorsF[plotrange],
        label = "ML")
    xlabel!("τ");ylabel!("gₐ");title!("Axial Charge")
    plot!(twinx(), xmirror=:true,grid=:false,ylims=(ytickvals[1],ytickvals[end]),
        xlims=(xtickvals[1],xtickvals[end]),xticks = (xtickvals,xtick0),
        yticks=(ytickvals,ytick0))
    annotate!(.34, 1.93, text("Axial Charge Estimates:",6, :left, :black))
    annotate!(.17, 1.91, text("Real: $sigChargeR($sigChargeSER) | χ²ᵥ = $chisqF ",6, :left, :purple))
    annotate!(.18, 1.89, text("ML: $sigChargeF($sigChargeSEF) | χ²ᵥ = $chisqR ",6, :left, :green))

    savefig("Axial Charge C3 Plot.png")

    return("Done with ML Axial!")
end
