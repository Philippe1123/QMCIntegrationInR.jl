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

    path=string(@__DIR__,"/res/Mapping/09092021/a=100/")
    for id = 1:18

        Data = load(string(path,id, "_Inv.jld2"))["data"]
        plotter(Data,path)

    end


end


function plotter(Data::Dict,path::String)

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


    figure()
    loglog(nbOfSamples, trueErrors, "-ks")
    #loglog(nbOfSamples, abs.(estimatedCubatureErrors), "-r*")
    loglog(nbOfSamples, abs.(trueCubatureErrors), "-rs")
    #loglog(nbOfSamples, abs.(estimatedTruncationErrors), "-g*")
    loglog(nbOfSamples, abs.(trueTruncationErrors), "-gs")
    loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
    loglog(nbOfSamples, nbOfSamples .^ -2, "--r")
    loglog(nbOfSamples, nbOfSamples .^ -3, "--y")
    grid(which="both",ls="-")


    str = string(
        QMCType,
        "\n",
        "s = ",
        s,
        ", shift = ",
        M,
        ", params = ",
        Data[13],
        ", alpha = ",
        Data[14],
    )

    strname = string(
        QMCType,
        "_",
        "s=",
        s,
        "_shift=",
        M,
        "_params=",
        Data[13],
        "_alpha=",
        Data[14],
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
