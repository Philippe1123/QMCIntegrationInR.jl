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
        Data = load(string(path,17, "_v12_res_tau=1.0_lat.jld2"))["data"]
       plotter(Data,path,"s",(255, 89, 72)./255)

        Data = load(string(path,17, "_v12_res_tau=1.0_lat_tent.jld2"))["data"]
       plotter(Data,path,"D",(255, 122, 40)./255)

        Data = load(string(path,18, "_v12_res_tau=1.0_sob.jld2"))["data"]
       plotter(Data,path,"^",(67, 133, 255)./255)

        Data = load(string(path,18, "_v12_res_tau=1.0_sob2.jld2"))["data"]
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
    ),fontsize=12)
    str = string(
    "s = ",
    s,
    ", shift = ",
    M,
    ", τ = ",
    "1.0"
)
    sz = 15
    title(str,fontsize=sz)

  #  end

  
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
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", τ = ",
"1.0"
)
sz = 15
title(str,fontsize=sz)
ax = gca()

ax[:yaxis][:set_ticks_position]("left")
ax[:yaxis][:set_label_position]("left")
plt.yticks(fontsize=sz)

PyPlot.savefig("tau1_PDE.pdf",transparent = "true",bbox_inches="tight")
##############################3

figure()
Data = load(string(path,17, "_v12_res_tau=2.0_lat.jld2"))["data"]
plotter(Data,path,"s",(255, 89, 72)./255)

Data = load(string(path,17, "_v12_res_tau=2.0_lat_tent.jld2"))["data"]
plotter(Data,path,"D",(255, 122, 40)./255)

Data = load(string(path,18, "_v12_res_tau=2.0_sob.jld2"))["data"]
plotter(Data,path,"^",(67, 133, 255)./255)

Data = load(string(path,18, "_v12_res_tau=2.0_sob2.jld2"))["data"]
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
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", τ = ",
"2.0"
)
sz = 15
title(str,fontsize=sz)

#  end


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
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", τ = ",
"2.0"
)
sz = 15
title(str,fontsize=sz)
ax = gca()

ax[:yaxis][:set_ticks_position]("right")
ax[:yaxis][:set_label_position]("right")
plt.yticks(fontsize=sz)

PyPlot.savefig("tau2_PDE.pdf",transparent = "true",bbox_inches="tight")

##############

figure()
Data = load(string(path,17, "_v12_res_tau=3.0_lat.jld2"))["data"]
plotter(Data,path,"s",(255, 89, 72)./255)

Data = load(string(path,17, "_v12_res_tau=3.0_lat_tent.jld2"))["data"]
plotter(Data,path,"D",(255, 122, 40)./255)

Data = load(string(path,18, "_v12_res_tau=3.0_sob.jld2"))["data"]
plotter(Data,path,"^",(67, 133, 255)./255)

Data = load(string(path,18, "_v12_res_tau=3.0_sob2.jld2"))["data"]
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
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", τ = ",
"3.0"
)
sz = 15
title(str,fontsize=sz)

#  end


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
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", τ = ",
"3.0"
)
sz = 15
title(str,fontsize=sz)
ax = gca()

ax[:yaxis][:set_ticks_position]("left")
ax[:yaxis][:set_label_position]("left")
plt.yticks(fontsize=sz)

PyPlot.savefig("tau3_PDE.pdf",transparent = "true",bbox_inches="tight")


##############

figure()
Data = load(string(path,17, "_v12_res_tau=4.0_lat.jld2"))["data"]
plotter(Data,path,"s",(255, 89, 72)./255)

Data = load(string(path,17, "_v12_res_tau=4.0_lat_tent.jld2"))["data"]
plotter(Data,path,"D",(255, 122, 40)./255)

Data = load(string(path,18, "_v12_res_tau=4.0_sob.jld2"))["data"]
plotter(Data,path,"^",(67, 133, 255)./255)

Data = load(string(path,18, "_v12_res_tau=4.0_sob2.jld2"))["data"]
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
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", τ = ",
"4.0"
)
sz = 15
title(str,fontsize=sz)

#  end


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
),fontsize=12)
str = string(
"s = ",
s,
", shift = ",
M,
", τ = ",
"4.0"
)
sz = 15
title(str,fontsize=sz)
ax = gca()

ax[:yaxis][:set_ticks_position]("right")
ax[:yaxis][:set_label_position]("right")
plt.yticks(fontsize=sz)

PyPlot.savefig("tau4_PDE.pdf",transparent = "true",bbox_inches="tight")

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
