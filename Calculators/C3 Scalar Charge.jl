using Plots
using Statistics
using LsqFit
include("C:\\Users\\Drew\\github\\SULI-LQCD\\myfunctions.jl")

Zsv = .877
ZsvSE = .001

cd("C:\\Users\\Drew\\github\\SULI-LQCD\\Data")
global datafile=open("C2ProtonGGData.txt","r");
sourcedata=readlines(datafile)
close(datafile)
sourcedata = split(split(split(sourcedata[5],"[")[2],"]")[1], ",")
sourcedata = [parse(Float64,i) for i in sourcedata]

#Starting the process at T0 folder in AMA
Tindex=0
dir0="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
cd(dir0)

if (isfile("C:\\Users\\Drew\\github\\SULI-LQCD\\Data\\C3ScalarChargeData.txt"))
    rm("C:\\Users\\Drew\\github\\SULI-LQCD\\Data\\C3ScalarChargeData.txt")
end

# gauge configuration list
gaugeconfigs = ["$(748 + 16*i)" for i in range(0,Integer((1420-748)/16),step=1)]
filter!(e->e∉["956","1004","1036","1052"], gaugeconfigs)
fpconfig = 112 #files per gauge config (112 measurements per gauge config)
# fpconfig is used to average all the gauge configs

# For loop Iterating over the "TX" files indicating source time
numconfigs = 39
timelength = 64
tempdata = []
binneddata = zeros(Float64,(numconfigs,timelength))
global C2 = []
@progress (name="Plotting...") for i in range(0,48,step=8)
    Tindex=i

    # directory for source-time folder
    dir="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
    cd(dir)
    # looping over all position (x,y,z) folders within TX folder
    for foldername in readdir()
        dir = string(dir, "\\", foldername)
        cd(dir)

        global filevector = zeros(Float64,(numconfigs,timelength))

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
            readuntil(test_file, "G0", keep = false)
            reldata = readuntil(test_file, "END_NUC3PT", keep=false)
            lines=split(reldata,"\n")

            # This is where finding C₂(t) begins #
            linesvector=[] # vector full of lines in the relevant data file

            # pushing data-containing split lines into a vector
            for l in range(6,length(lines)-5,step=1)
                push!(linesvector,split(lines[l]))
            end

            Op = zeros(Float64,timelength)

            # Fill operator vector with each line's relevant value (U-D)
            for l in range(1,length(linesvector),step=1)
                Op[l]=value(linesvector[l][2]) - value(linesvector[l][4])
            end
            Op = reshapevector(Op, Tindex) #Reshaping vectors so they are not time-shifted
            binneddata[rownum,:] += Op/fpconfig

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

# Turn binneddata into binned Jack replicates,
for i in range(1,timelength,step=1)
    binneddata[:,i] = Jackrep(binneddata[:,i])
end

#  Renormalizing by sourcedata
for i in range(1,numconfigs,step=1)
    binneddata[i,:] = binneddata[i,:]/(abs(sourcedata[i])*3.2)
end

cd("C:\\Users\\Drew\\github\\SULI-LQCD\\Data")
global datafile=open("gVRatioData.txt","r");
gV=readlines(datafile)
close(datafile)
gV = [split(split(split(gV[i],"[")[2],"]")[1], ",") for i in range(1,8,step=1)]

Ratio = zeros((length(gV[1]),length(gV)))
for i in range(1,length(Ratio[1,:]),step=1) # 1-8 tau
    for j in range(1,length(Ratio[:,1]),step=1) # 1-39 bins
        Ratio[j,i] = binneddata[j,i+1]/parse(Float64,gV[i][j])
    end
end
renormgs = [mean(Ratio[:,i]) for i in range(1,length(Ratio[1,:]),step=1)] # These are jack estimators
renormgsSE = [JackSE(Ratio[:,i]) for i in range(1,length(Ratio[1,:]),step=1)] # These are jack errors

physvalsSE = [sqrt( (renormgs[i]*ZsvSE)^2 + (Zsv*renormgsSE[i])^2 ) for i in range(1,length(renormgs),step=1)]
physvals = [Zsv * renormgs[i] for i in range(1,length(renormgs),step=1)]

Chargesren = [mean(Ratio[i,:]) for i in range(1,length(Ratio[:,1]),step=1)] # Jack replicates of charge
Chargeren = mean(Chargesren) #Jack estimator of Charge
ChargeSEren = JackSE(Chargesren)
chisqren = 0
for i in range(1,length(plateau),step=1)
    global chisqren += ((Chargeren-renormgs[i]))^2  / (renormgsSE[i]) / (length(renormgs)-1)
end
println("χ²/dof renorm fit = $chisqren")

# Populating Jack estimators and standard errors
stderrors = zeros(timelength)
finalvals = zeros(timelength)
for i in range(1,timelength,step=1)
    finalvals[i]=mean(binneddata[:,i])
    stderrors[i] = JackSE(binneddata[:,i])
end

plateau = 2:9
Charges = [mean(binneddata[:,i]) for i in plateau] # Jack replicates of charge
Charge = mean(Charges) #Jack estimator of Charge
ChargeSE = JackSE(Charges)
chisq = 0
for i in range(1,length(plateau),step=1)
    global chisq += ((Charge-finalvals[plateau[i]]))^2  / (stderrors[plateau[i]]) / (length(plateau)-1)
end
println("χ²/dof = $chisq")


cd("C:\\Users\\Drew\\github\\SULI-LQCD\\Data")
global dataoutfile = open("C3ScalarChargeData.txt","a") #saving data to file
write(dataoutfile,string(finalvals,"\n", stderrors, "\n", physvals, "\n", physvalsSE,
    "\n", chisq, "\n", Charge, "\n", ChargeSE,
    "\n", chisqren , "\n", Chargeren, "\n", ChargeSEren, "\n"))
close(dataoutfile)
