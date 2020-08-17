using Statistics
# Parses scientific value strings as floats
function value(str)
    parse(Float64,str)
end

# Finds the SE of a vector of Jackknife Expectation values (AKA Jackknife Replicates)
# Note: takes real parts of means and input values
function JackSE(JackrepVector)
    coefficient = ((length(JackrepVector)-1)/length(JackrepVector))
    JackSEsum=0
    for i in range(1,length(JackrepVector),step=1)
        JackSEsum = JackSEsum + (real(JackrepVector[i])-real(mean(JackrepVector)))^2
    end
    SE = sqrt(coefficient*JackSEsum)
    return(SE)
end

# Takes a vector of data and return a vector of Jackknife Replicates
function Jackrep(datavector)
    expectvector = []
    for i in range(1,length(datavector),step=1)
        expect = 0
        for j in range(1,length(datavector),step=1)
            if (j==i)
                continue
            else
                expect = expect + datavector[j]
            end
        end
        expect = expect / (length(datavector)-1)
        push!(expectvector,expect)
    end
    return(expectvector)
end

# Reshapes vectors how I want for zoomed plots
function reshapevector(vector, i)
    newvector = zeros(length(vector))
    for j in range(1,length(vector),step=1)
        if ((i+j) <= length(vector))
            newvector[j] = vector[i+j]
        else
            newvector[j] = vector[(i+j)%length(vector)]
        end
    end
    return (newvector)
end

# Creates m* vectors from C2 vectors, deletes outliers
function Emass(C2vector)
    Emassvector = []
    Emass = 0
    for j in range(1,length(C2vector),step=1)
        if (j==length(C2vector))
            Emass = C2vector[j]/C2vector[1]
        else
            Emass = C2vector[j]/C2vector[j+1]
        end
        push!(Emassvector,Emass)
    end
    Emassvector=real(log.(Complex.(Emassvector)))
    return(Emassvector)
end

# Converts lattice unit masses to MeV
function massconvert(mass,massSE)
    invspacing = 1.378
    invspacingSE = .007
    n_mass = mass*invspacing
    n_massSE = sqrt( (mass*invspacingSE)^2 + (invspacing*massSE)^2 )
    return([1000*n_mass,1000*n_massSE])
end

# Takes in matrices of axial charge replicates and vector charge replicates,
# produces matrix of physical axial charges and standard error vector.
function physax(Ax,Ve)
    PAx = zeros((length(Ax[:,1]),length(Ax[1,:])))
    SE = zeros(length(PAx[1,:]))
    for i in range(1,length(Ax[:,1]),step=1) #Iterate over all rows
        for j in range(1,length(Ax[1,:]),step=1) #Iterate over all columns
            PAx[i,j] = Ax[i,j]/Ve[i,j]
        end
    end
    SE = [JackSE(PAx[:,i]) for i in range(1,length(PAx[1,:]),step=1)]
    return(PAx, SE)
end

# Takes in matrices of scalar charge replicates and vector charge replicates,
# produces matrix of physical scalar charges and standard error vector.
function physsc(Sc,Ve)
    Zsv = .878 # Zs/Zv
    ZsvSE = .001# Error of Zs/Zv

    PSc = zeros((length(Sc[:,1]),length(Sc[1,:])))
    SE = zeros(length(PSc[1,:]))
    for i in range(1,length(Sc[:,1]),step=1) #Iterate over all rows
        for j in range(1,length(Sc[1,:]),step=1) #Iterate over all columns
            PSc[i,j] = Sc[i,j]/Ve[i,j]
        end
    end

    # Initial error w/o propagation
    SE = [JackSE(PSc[:,i]) for i in range(1,length(PSc[1,:]),step=1)]
    # Propagating Error
    SE = [sqrt( (mean(PSc[:,i])*ZsvSE)^2 + (Zsv*SE[i])^2 ) for i in range(1,length(SE),step=1)]
    Psc = Psc.*Zsv
    return(PSc, SE)
end

# Plots the physical axial charge with inputs of
# real axial replicates & error, fake axial replicates & error,
# and plotrange
function plotpax(PAxR, SER, PAxF, SEF, plotrange)
    # Changing directory
    plotdir = "C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\FinalPlots"
    cd(plotdir)
    # Setting plot boundaries
    xtickvals = [i for i in range(0,9,step=1)]
    ytickvals = [(1 +.05*i) for i in range(0,8,step=1)]
    ytick0 = ["" for i in range(1,length(ytickvals),step=1)]
    xtick0 = ["" for i in range(1,length(xtickvals),step=1)]

    # Plotting ML and Real comparison
    scatter(plotrange[1]-1:plotrange[end]-1,PAxR,marker=(:x),markercolor=(:purple),
        linecolor=(:purple),markerstrokecolor=(:purple),yerror=SER,
        dpi=600,grid=false,xlims=(xtickvals[1],xtickvals[end]),
        ylims=(ytickvals[1],ytickvals[end]),xticks=xtickvals,yticks=ytickvals,frame=(:box),
        foreground_color_legend = nothing, background_color_legend=nothing, label = "REAL",
        legend = ((.85,.95)), legendfontsize = 7)
    scatter!(plotrange[1]-1:plotrange[end]-1,PAxF,marker=(:x),markercolor=(:green),
        linecolor=(:green),markerstrokecolor=(:green),yerror=SEF,
        label = "ML")
    hline!(([1.2755]),linecolor=(:blue),label="")
    annotate!(0.05, 1.285, text("Experiment: 1.2732(23)",6, :left))
    hline!(([1.2709]),linecolor=(:blue),label="")
    xlabel!("τ");ylabel!("gₐZₐ");title!("Physical Axial Charge")
    savefig("Physical Axial Charge C3 Plot.png")
end

# Plots the physical scalar charge with inputs of
# real scalar replicates & error, fake scalar replicates & error,
# and plotrange
function plotpsc(PScR, SER, PScF, SEF, plotrange)
    # Changing directory
    plotdir = "C:\\Users\\Drew\\github\\SULI-LQCD\\ML\\FinalPlots"
    cd(plotdir)
    # Setting plot boundaries
    xtickvals = [i for i in range(0,9,step=1)]
    ytickvals = [(.2 +.2*i) for i in range(0,6,step=1)]
    ytick0 = ["" for i in range(1,length(ytickvals),step=1)]
    xtick0 = ["" for i in range(1,length(xtickvals),step=1)]

    # Plotting ML and Real comparison
    scatter(plotrange[1]-1:plotrange[end]-1,PScR,marker=(:x),markercolor=(:purple),
        linecolor=(:purple),markerstrokecolor=(:purple),yerror=SER,
        dpi=600,grid=false,xlims=(xtickvals[1],xtickvals[end]),
        ylims=(ytickvals[1],ytickvals[end]),xticks=xtickvals,yticks=ytickvals,frame=(:box),
        foreground_color_legend = nothing, background_color_legend=nothing, label = "REAL",
        legend = ((.85,.95)), legendfontsize = 7)
    scatter!(plotrange[1]-1:plotrange[end]-1,PScF,marker=(:x),markercolor=(:green),
        linecolor=(:green),markerstrokecolor=(:green),yerror=SEF,
        label = "ML")
    xlabel!("τ");ylabel!("gₛZₛ");title!("Physical Scalar Charge")
    savefig("Physical Scalar Charge C3 Plot.png")
end
