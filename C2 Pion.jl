using Plots
using Statistics
using LsqFit
#Starting the process at T0 folder in AMA
Tindex=0
dir0="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
cd(dir0)

if (isfile("C:\\Users\\Drew\\github\\SULI-LQCD\\C2PionData.txt"))
    rm("C:\\Users\\Drew\\github\\SULI-LQCD\\C2PionData.txt")
end

# Defining a function to go from scientific notation to floats
function value(str)
    parse(Float64,str)
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
        if ((i+j) <= length(vector))
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
    Emassvector=real(log.(Complex.(Emassvector)))
    return(Emassvector)
end

function EmassSE(C2vector,C2SEvector)
    EmassSEvector = []
    EmassSE = 0
    for i in range(1,length(C2vector),step=1)
        if (i==length(C2vector))
            a = (C2SEvector[i]/(C2vector[1]))
            b = (C2vector[i]*C2SEvector[1]/(C2vector[1]^2))
            EmassSE = sqrt(a^2 + b^2)*C2vector[1]/C2vector[i]
        else
            a = (C2SEvector[i]/(C2vector[i+1]))
            b = (C2vector[i]*C2SEvector[i+1]/(C2vector[i+1]^2))
            EmassSE = sqrt(a^2 + b^2)*C2vector[i+1]/C2vector[i]
        end
        push!(EmassSEvector,EmassSE)
    end
    EmassSEvector=real(log.(Complex.(EmassSEvector)))
    return(EmassSEvector)
end

################################################################################
########################### END OF FUNCTIONS ###################################
################################################################################

# For loop Iterating over the "TX" files indicating source time
numfiles=7*16*39
tempdata = []
alldata = zeros(ComplexF64,(numfiles,64))
global C2 = []
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

        global filevector = []

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
            readuntil(test_file, "GAM_5", keep = false)
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

            Op1 = zeros(ComplexF64,(length(times)))

            # Fill operator arrays
            for l in range(1,length(linematrices),step=1)
                Op1[l]=value(linematrices[l][2])+(value(linematrices[l][3]))im

            end

            push!(filevector,Op1[:])

        end

        # Rearrange Operatormatrix into C₂(μₐ,t)
        for k in range(1,length(filevector),step=1)
            gaugeconfig=[] # Slice t=l across all μₐ
            for l in range(1,length(filevector[1]),step=1)
                push!(gaugeconfig,filevector[k][l])
            end
            push!(C2,gaugeconfig)
        end

        for i in range(1,length(C2),step=1)
            push!(tempdata,C2[i])
        end

        # This is where finding C₂(t) ends #

        # Printing confirmation and path of completed file
        println("Done with $dir")

        # Resetting directory back to TX so that we can append x,y,z onto it
        # (Also so we don't get the wrong result from "readdir()" when we start again)
        dir="C:\\Users\\Drew\\Desktop\\BNL DATA\\AMA\\T$Tindex"
        cd(dir)
    end
end
close(test_file)
for i in range(1,length(C2),step=1)
    alldata[i,:]=C2[i]
end
alldata=real.(alldata)

# Binning data
binnedmeans=zeros((39,length(alldata[1,:])))
for i in range(1,39,step=1) # Iterate over gauge configurations
    binnedmatrix = zeros((7*16,64))
    for j in range(1,7*16,step=1)
        index = i+39(j-1)
        binnedmatrix[j,:]=alldata[index,:]
    end
    #Do jackknife error on binned data
    meanvector = zeros(length(binnedmatrix[1,:]))
    for j in range(1,length(binnedmatrix[1,:]),step=1)
        meanvector[j]=mean(binnedmatrix[:,j])
    end
    binnedmeans[i,:] = meanvector
end

# Turn binnedmeans into binned Jack replicates
for i in range(1,length(binnedmeans[1,:]),step=1)
    binnedmeans[:,i] = Jackrep(binnedmeans[:,i])
end

# Fitting data to find m*
fitmassreps = zeros((2,length(binnedmeans[:,1])))
plateau = 5:10
endplateau = 55:60
model(t,p) = p[1]*exp.(-p[2]*t)
covar = zeros((2,2))
for i in range(1,length(binnedmeans[:,1]),step=1) # folding data vvvv
    global fit = curve_fit(model,plateau,(binnedmeans[i,plateau]+reverse(binnedmeans[i,endplateau]))/2,[1,.1])
    global covar += estimate_covar(fit)
    fitmassreps[:,i] = fit.param
end
covar = covar / length(binnedmeans[:,1])
EffectiveMass = mean(fitmassreps[2,:])
Amplitude = mean(fitmassreps[1,:])
EffectiveMassSE = JackSE(fitmassreps[2,:])

Fitfunction(t) = Amplitude*ℯ^(-EffectiveMass*t)

# Populating final Jack estimators and errors
stderrors = zeros(length(binnedmeans[1,:]))
finalvals = zeros(length(binnedmeans[1,:]))

for i in range(1,length(binnedmeans[1,:]),step=1)
    stderrors[i] = JackSE(binnedmeans[:,i])
    finalvals[i]=mean(binnedmeans[:,i])
end

# Populating Jack replicates of effective mass
Effmassrep=zeros((39,length(alldata[1,:])))

for i in range(1,length(binnedmeans[:,1]),step=1)
    Effmassrep[i,:]=Emass(binnedmeans[i,:])
end

# Populating effective mass and error for m* plot
Effmass = zeros(length(Effmassrep[1,:]))
EffmassSE = zeros(length(Effmassrep[1,:]))

for i in range(1,length(Effmassrep[1,:]),step=1)
    Effmass[i]=mean(Effmassrep[:,i])
    EffmassSE[i] = JackSE(Effmassrep[:,i])
end

# Finding χ² of our fit
chisq = 0
for i in plateau
    global chisq +=  (((finalvals[i]+finalvals[64-i+1])/2 - Fitfunction(i))^2)/(((stderrors[i]+stderrors[64-i+1])/2)^2)
end
chisq = chisq / 2
println("χ²/dof = $chisq")

cd("C:\\Users\\Drew\\github\\SULI-LQCD")
global dataoutfile = open("C2PionData.txt","a") #saving data to file -> C2, C2 error, m* plot, m* error plot, m* estimate, m* estimate error
write(dataoutfile,string(finalvals,"\n", stderrors,
    "\n", Effmass, "\n", EffmassSE,
    "\n", EffectiveMass, "\n", EffectiveMassSE,
    "\n", "χ²=$chisq", "\n", "Covariance matrix: ", covar, "\n"))
close(dataoutfile)
