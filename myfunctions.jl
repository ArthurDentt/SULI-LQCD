using Statistics
# Parses scientific value strings as floats
function value(str)
    parse(Float64,str)
end

# Finds the SE of a vector of Jackknife Expectation values (AKA Jackknife Replicates)
# Note: takes real parts of means and values input
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

function vcovar(v1,v2)
    sum = 0
    for i in range(1,length(v1),step=1)
        sum += (v1[i]-mean(v1))*(v2[i]-mean(v2))/(length(v1)-1)
    end
    return(sum)
end

function mcovar(M)
    Mn = zeros((length(M[1,:]),length(M[1,:])))
    for i in range(1,length(M[1,:]),step=1) # iterate over columns
        for j in range(1,length(M[1,:]),step=1) # iterate over columns again
            Mn[i,j] = vcovar(M[:,i],M[:,j]) # covariance of columns i,j
        end
    end
    return(Mn)
end
