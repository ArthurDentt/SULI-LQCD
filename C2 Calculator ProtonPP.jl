using Plots
using Statistics

#Starting the process at T0 folder in AMA
Tindex=0
dir0="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
cd(dir0)

# Defining a function to go from a complex number z to |z|
# Maybe there is already one for this but I didn't find it lol
function modulus(z)
    sqrt(real(z)^2 + imag(z)^2)
end

# Defining a function to go from scientific notation to floats
function value(str)
    parts = split(str,"e")
    digits = parse(Float64,parts[1])
    power = parse(Int,parts[2])
    number=digits*(10.0^(power))
    return(number)
end

# Defining a function to find the SE of a vector of Jackknife Expectation values
# (AKA Jackknife Replicates) Note: takes real parts of means and values input
function JackSE(JackrepVector)
    coefficient = ((length(JackrepVector)-1)/length(JackrepVector))
    JackSEsum=0
    for i in range(1,length(JackrepVector),step=1)
        JackSEsum = JackSEsum + (real(JackrepVector[i])-real(mean(JackrepVector)))^2
    end
    SE = sqrt(coefficient*JackSEsum)
    return(SE)
end

# Defining a function to to take a vector of data and return a vector of
# Jackknife Replicates
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

# Defining a function to reshape vectors how I want for zoomed plots
function reshapevector(vector, i)
    newvector = zeros(length(vector))
    for j in range(1,length(vector),step=1)
        if (i+j <= length(vector))
            newvector[j] = vector[i+j]
        else
            newvector[j] = vector[(i+j)%length(vector)]
        end
    end
    return (newvector)
end

# Defining a function to create m* vectors from C2 vectors, deletes outliers
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
    return(Emassvector)
end

slices = [27:40,27:40,27:40,23:40,25:38,25:35,24:42]

# For loop Iterating over the "TX" files indicating source time
@progress (name="Plotting...") for i in range(0,48,step=8)
    Tindex=i
    # I'm going to leave my directory path here so you can see how this works
    # and what to change should be fairly clear
    dir="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
    cd(dir)
    # looping over all position (x,y,z) folders within TX folder
    for j in range(1,length(readdir()),step=1)
        rdir=readdir()[j]
        dir = string(dir, "\\", rdir)
        cd(dir)

        Operatormatrix = [] # It's a matrix full of operator values
        C2 = [] # Defining C₂(t,μₐ) -> μₐ is Gauge Config α
        realEC2=[] # Real part of expectation value of C₂(t)
        jackreplicates = [] # matrix of jackknife replicates
        stderror = [] # error vector sqrt(variance(jackknife))
        global jackestimates = []

        # Starting the main loop running over each text file ending in .dat.###
        for k in range(0,42,step=1)
            findex = 748 + k * 16
            try
                global test_file=open("nuc3pt.dat.$findex","r+");
            catch (SystemError)
                # println("No file named \"nuc3pt.dat.$findex\" exists.")
                # Un-comment to see what files are missing, I don't really care though
                continue
            end

            # Reading to keyword, splitting lines, instantiating + zeroing everything
            readuntil(test_file, "NUC_G5C_NP", keep = false)
            reldata = readuntil(test_file, "ENDPROP", keep=false)
            lines=split(reldata,"\n")

            # This is where finding C₂(t) begins #
            times=[]; # times
            linematrices=[] # matrix full of pieces of each line in the relevant data file

            # pushing split lines into a matrix
            for l in range(2,length(lines)-1,step=1)
                push!(linematrices,split(lines[l]))
                push!(times,parse(Int,linematrices[l-1][1]))
            end

            # Define operator arrays
            Op1 = zeros(ComplexF64,(length(times)))

            # For now I am not using Operator 2 (it is the negative parity C₂(t))
            # Un-comment all similar lines to use it
            #Op2 = zeros(ComplexF64,(length(times)))

            # Fill operator arrays
            for l in range(1,length(linematrices),step=1)
                Op1[l]=value(linematrices[l][2])+(value(linematrices[l][3]))im

                # For now I am not using Operator 2
                #Op2[i]=value(linematrices[i][4])+(value(linematrices[i][5]))im
            end

            push!(Operatormatrix,Op1[:])
            #println("Pushed Op1 to Operatormatrix, Filename: \"nuc3pt.dat.$findex\"")

        end

        # Rearrange Operatormatrix into C₂(t,μₐ)
        for k in range(1,length(Operatormatrix[1]),step=1)
            timeslice=[] # Slice t=l across all μₐ
            for l in range(1,length(Operatormatrix),step=1)
                push!(timeslice,Operatormatrix[l][k])
            end
            push!(C2,timeslice)
        end

        # Filling jackreplicates matrix
        for k in range(1,length(C2),step=1)
            jackknifevector=Jackrep(C2[k])
            push!(jackreplicates,jackknifevector)
        end

        # Filling Standard error matrix
        for k in range(1,length(jackreplicates),step=1)
            SE = JackSE(jackreplicates[k])
            push!(stderror,SE)
        end

        # Filling jackknife estimate vector
        for k in range(1,length(jackreplicates),step=1)
            push!(jackestimates,real(mean(jackreplicates[k])))
        end

        if (j==1)

            slice = 19:length(jackestimates)-19
            resjack = reshapevector(jackestimates,i)
            reserror = reshapevector(stderror,i)

            #             zoomed-in plot                #
            #zoomplot=plot(slice,resjack[slice],marker=(:circle),legend=false,grid=false,yerror=reserror[slice])
            #xlabel!("t -> shifted left by $i");ylabel!("Re(<C₂>)");title!("Proton-PP Re(<C₂>)(t)");
            # Displaying plots takes a while, so I don't do it
            #display(mec2plot)
            #cd("C:\\Users\\Drew\\Desktop\\C2Graphs\\Proton_PP\\RealEC2(t)")
            #savefig("Zoomed in - T=$(Tindex), $(rdir).png")
            #global dataoutfile = open("C2ProtonOutput.txt","a")
            #write(dataoutfile,string(resjack[slices[Int((i/8)+1)]],"\n"))   # Not Used for now #
            #close(dataoutfile)

            #             Ratio plot C2[t] / C2[t+1]                            #
            #Emassplot=plot(1:length(Emass(jackestimates)),Emass(jackestimates),marker=(:circle),legend=false)
            #xlabel!("t");ylabel!("m*");title!("Proton PP m*(t)")
            # Displaying plots takes a while, so I don't do it
            #display(mec2plot)
            #cd("C:\\Users\\Drew\\Desktop\\C2Graphs\\Proton_PP\\RealEC2(t)")
            #savefig("EFFECTIVE MASS T=$(Tindex), $(rdir).png")
            #println("Emass: $(Emass(jackestimates))")
        end


        #             Plot Real Part With Error                 #
        rec2plot=plot(1:length(jackestimates),jackestimates,marker=(:circle),legend=false,grid=false,yerror=stderror)
        xlabel!("t");ylabel!("Re(<C₂>)");title!("Proton-PP Re(<C₂>)(t)");
        # Displaying plots takes a while, so I don't do it
        #display(mec2plot)
        cd("C:\\Users\\Drew\\Desktop\\C2Graphs\\Proton_PP\\RealEC2(t)")
        savefig("T=$(Tindex), $(rdir).png")


        # This is where finding C₂(t) ends #

        # Printing confirmation and path of completed file
        println("Done with $dir")

        # Resetting directory back to TX so that we can append x,y,z onto it
        # (Also so we don't get the wrong result from "readdir()" when we start again)
        dir="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
        cd(dir)
    end
end

# Note: because this is all in a big loop, you won't see anything in the workspace
# variable explorer :(
