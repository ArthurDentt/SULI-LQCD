dir0="C:\\Users\\Drew\\Desktop\\FakeData"
cd(dir0)
filelist = readdir()
for i in filelist
    println(i)
    if (!(occursin("t",i[9:end-3])))
        rm("C:\\Users\\Drew\\Desktop\\FakeData\\$i")
    end
end
