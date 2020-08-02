using Plots
using Statistics
using LsqFit
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")
#Starting the process at T0 folder in AMA
Tindex=0
dir0="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
cd(dir0)

if (isfile("C:\\Users\\Drew\\github\\SULI-LQCD\\Data\\C2PionData.txt"))
    rm("C:\\Users\\Drew\\github\\SULI-LQCD\\Data\\C2PionData.txt")
end

# data-saving directory
datadir = "C:\\Users\\Drew\\github\\SULI-LQCD\\Data"

# gauge configuration list
gaugeconfigs = ["$(748 + 16*i)" for i in range(0,Integer((1420-748)/16),step=1)]
filter!(e->e∉["956","1004","1036","1052"], gaugeconfigs)
fpconfig = 112

# For loop Iterating over the "TX" files indicating source time
numconfigs = 39
timelength = 64
tempdata = []
binneddata = zeros(Float64,(numconfigs,64))
global C2 = []
@progress (name="Plotting...") for i in range(0,48,step=8)
    Tindex=i
    # I'm going to leave my directory path here so you can see how this works
    # and what to change should be fairly clear
    dir="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
    cd(dir)
    # looping over all position (x,y,z) folders within TX folder
    for foldername in readdir()
        dir = string(dir, "\\", foldername)
        cd(dir)

        global filevector = zeros(Float64,(numconfigs,64))

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

            #Finding row number (based on gauge config) 748->1, ..., 1420->39
            rownum = 0
            for gaugenum in gaugeconfigs
                rownum += 1
                if (occursin(gaugenum,"$findex"))
                    break
                end
            end

            # Reading to keyword, splitting lines, instantiating + zeroing everything
            readuntil(test_file, "GAM_5", keep = false)
            reldata = readuntil(test_file, "ENDPROP", keep=false)
            lines=split(reldata,"\n")

            # This is where finding C₂(t) begins #
            linesvector=[] # vector full of lines in the relevant data file

            # pushing data-containing split lines into a vector
            for l in range(2,length(lines)-1,step=1)
                push!(linesvector,split(lines[l]))
            end

            Op1 = zeros(Float64,(length(lines)-2))

            # Fill operator vector with each line's relevant value
            for l in range(1,length(linesvector),step=1)
                Op1[l]=value(linesvector[l][2])
            end

            #Op1 = reshapevector(Op1, Tindex) #Do we not need this???
            binneddata[rownum,:] += Op1/fpconfig

        end
        # Printing confirmation and path of completed file
        println("Done with $dir")

        # Resetting directory back to TX so that we can append x,y,z onto it
        # (Also so we don't get the wrong result from "readdir()" when we start again)
        dir="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
        cd(dir)
    end
end
close(test_file)

#folding binned data before it becomes jackknife replicates
foldsavedbins = (binneddata[:,:]+reverse(binneddata[:,:],dims=2))/2

# Turn binneddata into binned Jack replicates, Populating final Jack estimators and errors
stderrors = zeros(length(binneddata[1,:]))
finalvals = zeros(length(binneddata[1,:]))
for i in range(1,length(binneddata[1,:]),step=1)
    binneddata[:,i] = Jackrep(binneddata[:,i])
    finalvals[i]=mean(binneddata[:,i])
    stderrors[i] = JackSE(binneddata[:,i])
end

# Fitting data to find m*
fitmassreps = zeros((2,length(binneddata[:,1])))
plateau = 7:11
model(t,p) = p[1]*exp.(-p[2]*t)
foldbins = (binneddata[:,:]+reverse(binneddata[:,:],dims=2))/2
for i in range(1,length(binneddata[:,1]),step=1) # folding data vvvv
    global fit = curve_fit(model,plateau,foldbins[i,plateau],[1,.01])
    fitmassreps[:,i] = fit.param
end
EffectiveMass = mean(fitmassreps[2,:])
Amplitude = mean(fitmassreps[1,:])
EffectiveMassSE = JackSE(fitmassreps[2,:])

Fitfunction(t) = Amplitude*ℯ^(-EffectiveMass*t)

covariancemat = cov(foldsavedbins[:,plateau])
icov = inv(covariancemat)
foldfinalvals = (finalvals + reverse(finalvals))/2
# Finding χ² of our fit
chisq = 0
for i in range(1,length(plateau),step=1)
    for j in range(1,length(plateau),step=1)
        global chisq += (( (-foldfinalvals[plateau[i]] + Fitfunction(plateau[i])) )*
            icov[i,j]*( (-foldfinalvals[plateau[j]] + Fitfunction(plateau[j])) )) / (length(plateau)-2)
    end
end
println("χ²/dof = $chisq")

#Changing to data directory and ouputting useful info for plotting / calcs.
cd(datadir)
global dataoutfile = open("C2PionData.txt","a") #saving data to file -> C2, C2 error, m* plot, m* error plot, m* estimate, m* estimate error
write(dataoutfile,string(finalvals,"\n", stderrors,
    "\n", EffectiveMass, "\n", EffectiveMassSE,
    "\n", "χ²=$chisq", "\n"))
close(dataoutfile)
