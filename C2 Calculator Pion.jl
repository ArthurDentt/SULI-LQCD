using Plots
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

# For loop Iterating over the "TX" files indicating source time
for i in range(0,48,step=8)
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
        jackreplicates = []
        stderror = []
        jackestimates = []

        # Starting the main loop running over each text file ending in .dat.###
        for i in range(0,42,step=1)
            findex = 748 + i * 16
            try
                global test_file=open("nuc3pt.dat.$findex","r+");
            catch (SystemError)
                # println("No file named \"nuc3pt.dat.$findex\" exists.")
                # Un-comment to see what files are missing, I don't really care though
                continue
            end

            # Reading to keyword, splitting lines, instantiating + zeroing everything
            readuntil(test_file, "GAM_5", keep = false)
            reldata = readuntil(test_file, "ENDPROP", keep=false)
            lines=split(reldata,"\n")

            # This is where finding C₂(t) begins #
            times=[]; # times
            linematrices=[] # matrix full of pieces of each line in the relevant data file

            # pushing split lines into a matrix
            for i in range(2,length(lines)-1,step=1)
                push!(linematrices,split(lines[i]))
                push!(times,parse(Int,linematrices[i-1][1]))
            end

            # Define operator arrays
            Op1 = zeros(ComplexF64,(length(times)))

            # Fill operator arrays
            for i in range(1,length(linematrices),step=1)
                Op1[i]=value(linematrices[i][2])+(value(linematrices[i][3]))im

            end

            push!(Operatormatrix,Op1[:])
            #println("Pushed Op1 to Operatormatrix, Filename: \"nuc3pt.dat.$findex\"")

        end

        # Rearrange Operatormatrix into C₂(t,μₐ)
        for i in range(1,length(Operatormatrix[1]),step=1)
            timeslice=[] # Slice t=i across all μₐ
            for j in range(1,length(Operatormatrix),step=1)
                push!(timeslice,Operatormatrix[j][i])
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

        #             Plot Real Part                            #
        rec2plot=plot(1:length(jackestimates),jackestimates,marker=(:circle),legend=false,yerror=stderror)
        xlabel!("t");ylabel!("Re(<C₂>)");title!("Pion-GAM_5 Re(<C₂>)(t)")
        # Displaying plots takes a while, so I don't do it
        #display(mec2plot)
        cd("C:\\Users\\Drew\\Desktop\\C2Graphs\\Pion\\RealEC2(t)")
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
