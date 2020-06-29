using Plots
using Statistics
#Starting the process at T0 folder in AMA
Tindex=0
println("Beginning output...")
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

# Defining a function to check if a string starts with an Integer
function startswithint(str)
    parts=split(str)
        try
            parse(Int,parts[1])
            return(true)
        catch ArgumentError
            return(false)
        end
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
@progress (name= "Scanning through files..." ) for i in range(0,48,step=8)
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
        C3 = [] # Defining C₃(t,μₐ) -> μₐ is Gauge Config α
        jackreplicates = []
        jackestimates = []
        stderror = []

        # Starting the main loop running over each text file ending in .dat.###
        @progress (name = "Gauge config $i, Position $j, T=$Tindex...") for i in range(0,42,step=1)
            findex = 748 + i * 16
            try
                global test_file=open("nuc3pt.dat.$findex","r+");
            catch (SystemError)
                # println("No file named \"nuc3pt.dat.$findex\" exists.")
                # Un-comment to see what files are missing, I don't really care though
                continue
            end

            # Reading to keyword, splitting lines, instantiating + zeroing everything
            readuntil(test_file, "G0", keep = false)
            reldata = readuntil(test_file, "END_NUC3PT", keep=false)
            lines=split(reldata,"\n")
            # This is where finding C₃(t) begins #
            times=[]; # times
            linematrices=[] # matrix full of pieces of each line in the relevant data file

            # pushing split lines into a matrix
            for i in range(6,length(lines)-1,step=1)
                push!(linematrices,split(lines[i]))
                push!(times,parse(Int,linematrices[i-5][1]))
            end

            # Define operator arrays
            Op1 = zeros(ComplexF64,(length(times)))

            # For now I am not using Operator 2 (it is the negative parity C₃(t))
            # Un-comment all similar lines to use it
            #Op2 = zeros(ComplexF64,(length(times)))

            # Fill operator arrays
            for i in range(1,length(linematrices),step=1)
                Op1[i]=value(linematrices[i][2])+(value(linematrices[i][3]))im

                # For now I am not using Operator 2
                #Op2[i]=value(linematrices[i][4])+(value(linematrices[i][5]))im
            end

            push!(Operatormatrix,Op1[:])
            #println("Pushed Op1 to Operatormatrix, Filename: \"nuc3pt.dat.$findex\"")

        end

        # Rearrange Operatormatrix into C₃(t,μₐ)
        for i in range(1,length(Operatormatrix[1]),step=1)
            timeslice=[] # Slice t=i across all μₐ
            for j in range(1,length(Operatormatrix),step=1)
                push!(timeslice,Operatormatrix[j][i])
            end

            push!(C3,timeslice)
        end

        # Filling jackreplicates matrix
        for k in range(1,length(C3),step=1)
            jackknifevector=Jackrep(C3[k])
            push!(jackreplicates,jackknifevector)
        end

        # Filling Standard error matrix
        for k in range(1,length(jackreplicates),step=1)
            SE = JackSE(jackreplicates[k])
            push!(stderror,SE)
        end

        # Filling jackknife estimate vector
        for k in range(1,length(jackreplicates),step=1)
            push!(jackestimates,real(mean(jackreplicates[k]))/8768.7709394415)
        end

        #             Plot Real Part                            #
        reC3plot=plot(1:length(jackestimates),jackestimates,marker=(:circle),legend=false)
        xlabel!("t");ylabel!("Re(<C₃>)");title!("Scalar Charge Re(<C₃>)(t)")
        # Displaying plots takes a while, so I don't do it
        #display(meC3plot)
        cd("C:\\Users\\Drew\\Desktop\\C3Graphs\\Scalar Charge\\RealEC3(t)")
        savefig("T=$(Tindex), $(rdir).png")

        # This is where finding C₃(t) ends #

        # Printing confirmation and path of completed file
        println("Done with $dir")

        # Resetting directory back to TX so that we can append x,y,z onto it
        # (Also so we don't get the wrong result from "readdir()" when we start again)
        dir="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
        cd(dir)
    end
end

println("Done!")
#Testing adding text
# Note: because this is all in a big loop, you won't see anything in the workspace
# variable explorer :(
