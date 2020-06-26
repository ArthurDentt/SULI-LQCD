cd("C:\\Users\\Drew\\Desktop\\Julia Files\\Project work\\Text Processing\\")
test_file=open("nuc3pt.dat.748","r+");
startword=true;
numericaldata=[];
linematrices=[]
# Defining the value function to interpret scientific notation
function value(str)
    parts = split(str,"e")
    digits = parse(Float64,parts[1])
    power = parse(Int,parts[2])
    number=digits*(10.0^(power))
    return(number)
end
# Defining a function to check if a string starts with an integer
function startswithint(str)
    parts=split(str)
        try
            parse(Int,parts[1])
            return(true)
        catch ArgumentError
            return(false)
        end
end
println("------------ BEGIN OUTPUT ----------")
# Starting the main while loop which iterates while Julia has found
# The keyword "Start" indicating that there is input to read
while (startword==true)
    numericaldata=[];
    startword = false;
    startcheck = readuntil(test_file, "START", keep = true)
    println(startcheck)
    if (length(startcheck) != 0)
    startword = true
        for i in range(1,100;step=1)
            currentline=readline(test_file)
            # println(currentline)
            if (startswithint(currentline)==true)
                push!(numericaldata, currentline)
                # println("Pushed: $currentline")
            elseif (startswith(currentline,"END") == true)
                break
            else
                println(currentline)
            end
        end
        times=zeros(1)
        global linematrices=[]
        Op1expect=0;
        Op2expect=0;
        #
        for i in range(1,length(numericaldata),step=1)
            push!(linematrices,split(numericaldata[i]))
            push!(times,parse(Int,linematrices[i][1]))
        end
        # Define operator arrays
        Op1 = zeros(ComplexF64,(length(times)))
        Op2 = zeros(ComplexF64,(length(times)))
        Op1expect=0
        Op2expect=0
        # Fill operator arrays
        for i in range(1,length(linematrices),step=1)
            #print(linematrices[i])
            Op1[i]=value(linematrices[i][2])+(value(linematrices[i][3]))im
            Op1expect=Op1expect+Op1[i]
            if (length(linematrices[i]) == 5)
                Op2[i]=value(linematrices[i][4])+(value(linematrices[i][5]))im
                Op2expect=Op2expect+Op2[i]
            end
        end
        Op1expect=Op1expect/length(times)
        Op2expect=Op2expect/length(times)
        println("Operator 1 Time average: $Op1expect")
        if (length(linematrices[1]) == 5)
            println("Operator 2 Time average: $Op2expect")
        end
    else
        break
    end
end
