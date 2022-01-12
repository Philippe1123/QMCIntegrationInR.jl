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
 #   path=string(@__DIR__,"/")
    path=string(@__DIR__,"/res-NoCommit/OutFor36Months/")
 #   path=string(@__DIR__,"/res/PDE/06122021_11elem/")
 #   path=string(@__DIR__,"/res/PDE/07122021_95elem/")

  figure()
    for id in [6]

        Data = load(string(path,id, "_v13_22_res.jld2"))["data"]
        plotter(Data,path)

    end


end


function plotter(Data::Dict,path::String;type_name::String="Nil")

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
    
    loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
    loglog(nbOfSamples, nbOfSamples .^ -2, "--r")
    loglog(nbOfSamples, nbOfSamples .^ -3, "--y")

    ylim([10^(-17),1])
    
    loglog(nbOfSamples, trueErrors, "-v")

    grid(which="both",ls="-")


    str = string(
        QMCType," ",type_name,
        "\n",
        "s = ",
        s,
        ", shift = ",
        M,
        ", params = ",
        Data[13],
        ", alpha = ",
        Data[14],
        "_a=",
        Data[16],
    )

    strname = string(
        QMCType,"_",type_name,
        "_",
        "s=",
        s,
        "_shift=",
        M,
        "_params=",
        Data[13],
        "_alpha=",
        Data[14],
        "_a=",
        Data[16],
    )
    title(str)
    ylabel("Rel. error")
    xlabel("N")
    legend((
        "true error",
        "true cubature error",
        "true truncation error",
        "N^-1",
        "N^-2",
        "N^-3",
    ))

    PyPlot.savefig(string(path,strname, ".png"))
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
