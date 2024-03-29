using Statistics: mean, std
using SpecialFunctions: erf, erfinv, gamma, gamma_inc
using StringLiterals
using PrettyTables
using Random
using DelimitedFiles
using JLD2
using FileIO
using PyPlot
using DataFrames


function main()
    path=string(@__DIR__,"/")
 #   path=string(@__DIR__,"/res-NoCommit/OutFor36Months/")
 #   path=string(@__DIR__,"/res/PDE/06122021_11elem/")
 #   path=string(@__DIR__,"/res/PDE/07122021_95elem/")

  #figure()
  #  for id in [15]
  figure()
        Data = load(string(path,15, "_v13_res0.5.jld2"))["data"]
       plotter(Data,path,"s",(255, 89, 72)./255)

        Data = load(string(path,15, "_v13_res0.5_tent.jld2"))["data"]
       plotter(Data,path,"D",(255, 122, 40)./255)

        Data = load(string(path,16, "_v13_res0.5_sob.jld2"))["data"]
       plotter(Data,path,"^",(67, 133, 255)./255)

        Data = load(string(path,16, "_v13_res0.5_sob2.jld2"))["data"]
       plotter(Data,path,"o",(0, 179, 89)./255)
       nbOfSamples = Data[11]
       s = Data[9]
       M = Data[10]

    loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
    loglog(nbOfSamples, nbOfSamples .^ -2, "--r")

    
    legend((
        "rank-1 Lattice",
        "Tent transformed rank-1 lattice",
        "Sobol'",
        "Sobol' Factor 2 Interlacing",
        "N^-1",
        "N^-2",
        "QMC std",
    ),fontsize=12)
    str = string(
    "s = ",
    s,
    ", shift = ",
    M,
    ", domain = ",
    "[",-0.5,",",0.5,"]"
)
    sz = 15
    title(str,fontsize=sz)
    PyPlot.savefig("05.pdf",transparent = "true",bbox_inches="tight")

  #  end

  ###########################################################################################

  figure()
  Data = load(string(path,15, "_v13_res1.0.jld2"))["data"]
 plotter(Data,path,"s",(255, 89, 72)./255)

  Data = load(string(path,15, "_v13_res1.0_tent.jld2"))["data"]
 plotter(Data,path,"D",(255, 122, 40)./255)

  Data = load(string(path,16, "_v13_res1.0_sob.jld2"))["data"]
 plotter(Data,path,"^",(67, 133, 255)./255)

  Data = load(string(path,16, "_v13_res1.0_sob2.jld2"))["data"]
 plotter(Data,path,"o",(0, 179, 89)./255)
 nbOfSamples = Data[11]
 s = Data[9]
 M = Data[10]

loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
loglog(nbOfSamples, nbOfSamples .^ -2, "--r")


legend((
  "rank-1 Lattice",
  "Tent transformed rank-1 lattice",
  "Sobol'",
  "Sobol' Factor 2 Interlacing",
  "N^-1",
  "N^-2",
  "QMC std",
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", domain = ",
"[",-1.0,",",1.0,"]"
)
sz = 15
title(str,fontsize=sz)
ax = gca()

ax[:yaxis][:set_ticks_position]("right")
ax[:yaxis][:set_label_position]("right")
plt.yticks(fontsize=sz)

PyPlot.savefig("10.pdf",transparent = "true",bbox_inches="tight")



###############################################################################3

figure()
Data = load(string(path,15, "_v13_res2.0.jld2"))["data"]
plotter(Data,path,"s",(255, 89, 72)./255)

Data = load(string(path,15, "_v13_res2.0_tent.jld2"))["data"]
plotter(Data,path,"D",(255, 122, 40)./255)

Data = load(string(path,16, "_v13_res2.0_sob.jld2"))["data"]
plotter(Data,path,"^",(67, 133, 255)./255)

Data = load(string(path,16, "_v13_res2.0_sob2.jld2"))["data"]
plotter(Data,path,"o",(0, 179, 89)./255)
nbOfSamples = Data[11]
s = Data[9]
M = Data[10]

loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
loglog(nbOfSamples, nbOfSamples .^ -2, "--r")


legend((
"rank-1 Lattice",
"Tent transformed rank-1 lattice",
"Sobol'",
"Sobol' Factor 2 Interlacing",
"N^-1",
"N^-2",
"QMC std",
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", domain = ",
"[",-2.0,",",2.0,"]"
)
sz = 15
title(str,fontsize=sz)
PyPlot.savefig("20.pdf",transparent = "true",bbox_inches="tight")

#  end

###########################################################################################

figure()
Data = load(string(path,15, "_v13_res3.0.jld2"))["data"]
plotter(Data,path,"s",(255, 89, 72)./255)

Data = load(string(path,15, "_v13_res3.0_tent.jld2"))["data"]
plotter(Data,path,"D",(255, 122, 40)./255)

Data = load(string(path,16, "_v13_res3.0_sob.jld2"))["data"]
plotter(Data,path,"^",(67, 133, 255)./255)

Data = load(string(path,16, "_v13_res3.0_sob2.jld2"))["data"]
plotter(Data,path,"o",(0, 179, 89)./255)
nbOfSamples = Data[11]
s = Data[9]
M = Data[10]

loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
loglog(nbOfSamples, nbOfSamples .^ -2, "--r")


legend((
"rank-1 Lattice",
"Tent transformed rank-1 lattice",
"Sobol'",
"Sobol' Factor 2 Interlacing",
"N^-1",
"N^-2",
"QMC std",
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", domain = ",
"[",-3.0,",",3.0,"]"
)
sz = 15
title(str,fontsize=sz)
ax = gca()

ax[:yaxis][:set_ticks_position]("right")
ax[:yaxis][:set_label_position]("right")
plt.yticks(fontsize=sz)

PyPlot.savefig("30.pdf",transparent = "true",bbox_inches="tight")
###############################



figure()
Data = load(string(path,15, "_v13_res4.0.jld2"))["data"]
plotter(Data,path,"s",(255, 89, 72)./255)

Data = load(string(path,15, "_v13_res4.0_tent.jld2"))["data"]
plotter(Data,path,"D",(255, 122, 40)./255)

Data = load(string(path,16, "_v13_res4.0_sob.jld2"))["data"]
plotter(Data,path,"^",(67, 133, 255)./255)

Data = load(string(path,16, "_v13_res4.0_sob2.jld2"))["data"]
plotter(Data,path,"o",(0, 179, 89)./255)
nbOfSamples = Data[11]
s = Data[9]
M = Data[10]

loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
loglog(nbOfSamples, nbOfSamples .^ -2, "--r")


legend((
"rank-1 Lattice",
"Tent transformed rank-1 lattice",
"Sobol'",
"Sobol' Factor 2 Interlacing",
"N^-1",
"N^-2",
"QMC std",
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", domain = ",
"[",-4.0,",",4.0,"]"
)
sz = 15
title(str,fontsize=sz)
PyPlot.savefig("40.pdf",transparent = "true",bbox_inches="tight")

#  end

###########################################################################################

figure()
Data = load(string(path,15, "_v13_res5.0.jld2"))["data"]
plotter(Data,path,"s",(255, 89, 72)./255)

Data = load(string(path,15, "_v13_res5.0_tent.jld2"))["data"]
plotter(Data,path,"D",(255, 122, 40)./255)

Data = load(string(path,16, "_v13_res5.0_sob.jld2"))["data"]
plotter(Data,path,"^",(67, 133, 255)./255)

Data = load(string(path,16, "_v13_res5.0_sob2.jld2"))["data"]
plotter(Data,path,"o",(0, 179, 89)./255)
nbOfSamples = Data[11]
s = Data[9]
M = Data[10]

loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
loglog(nbOfSamples, nbOfSamples .^ -2, "--r")


legend((
"rank-1 Lattice",
"Tent transformed rank-1 lattice",
"Sobol'",
"Sobol' Factor 2 Interlacing",
"N^-1",
"N^-2",
"QMC std",
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", domain = ",
"[",-5.0,",",5.0,"]"
)
sz = 15
title(str,fontsize=sz)
ax = gca()

ax[:yaxis][:set_ticks_position]("right")
ax[:yaxis][:set_label_position]("right")
plt.yticks(fontsize=sz)

PyPlot.savefig("50.pdf",transparent = "true",bbox_inches="tight")

#  end


end


function plotter(Data::Dict,path::String,type::String,tp::Tuple;type_name::String="Sequence")

    timings = Data[2]
    trueErrors = Data[3]
    estimatedCubatureErrors = Data[4]
    trueCubatureErrors = Data[5]
    estimatedTruncationErrors = Data[6]
    trueTruncationErrors = Data[7]
    QMCType = Data[8]
    s = Data[9]
    M = Data[10]
    nbOfSamples = Data[11]


  
    
    #loglog(nbOfSamples, abs.(estimatedCubatureErrors), "-r*")
    #loglog(nbOfSamples, abs.(trueCubatureErrors), ":rs",markerfacecolor="None",markeredgecolor="red")
    #loglog(nbOfSamples, abs.(estimatedTruncationErrors), "-g*")
    #loglog(nbOfSamples, abs.(trueTruncationErrors), "-.gs")
    
    #loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
    #loglog(nbOfSamples, nbOfSamples .^ -2, "--r")
    #loglog(nbOfSamples, nbOfSamples .^ -3, "--y")

    ylim([10^(-17),1])
    
    loglog(nbOfSamples, trueErrors, color=tp,marker=type)

    grid(which="both",ls="-")

    sz = 15


    ylabel("error",fontsize=sz)
    xlabel("Number of samples N",fontsize=sz)
    """
    legend((
        "N^-1",
        "N^-2",
        "N^-3",
        "QMC std",
    ),fontsize=sz)
    """
    plt.xticks(fontsize=sz)
    plt.yticks(fontsize=sz)
    
    ax = gca()
    
    ax[:yaxis][:set_ticks_position]("left")
    ax[:yaxis][:set_label_position]("left")
    
    plt.yticks(fontsize=sz)

    #
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")








end




main()
