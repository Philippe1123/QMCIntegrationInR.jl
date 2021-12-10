
using StringLiterals
using PrettyTables
using Random
using DelimitedFiles
using JLD2
using FileIO
using StatsPlots
using DataFrames
using GR
using SpecialFunctions
using Statistics

analytical_sol(a::Real, s::Int64, sigma::Real) =
    (
        (
            gamma((1 + sigma) / 2) -
            gamma((1 + sigma) / 2) * gamma_inc((1 + sigma) / 2, a^2 / 2, 0)[2]
        ) * 2^(sigma / 2) / sqrt(pi) + erf(a / sqrt(2))
    )^s


function main()
    path=string(@__DIR__)
    ##for id = 1:18

        Data = load(string(path,"/","18_v51_res.jld2"))["data"]
        plotter(Data,path)

    ##end


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
        "/",QMCType,
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




    data11=Data[11]
    data15=Data[15]


    data11=log.(Data[11])./log(2)
    data11=repeat(data11',M,1)
    data11=reshape(data11,1,size(data11,1)*size(data11,2))

    data15=reshape(data15,size(data15,1)*size(data15,2),1)
    gr(size = (1000, 1000))
    df = DataFrame(sample = vec(data11), res= vec(data15))
    @df df violin(:sample,:res,linewidth=0,xticks=log.(Data[11])./log(2),xlabel="log2 of samples",ylabel="expected value",label="")
    @df df boxplot!(:sample,:res,linewidth=1.5,label="",fillalpha=0.75)
    @df df dotplot!(:sample,:res,markersize = 2.5,markerstrokewidth = 1.5,markerstrokealpha = 0.2, markerstrokecolor = :black,xticks=log.(Data[11])./log(2),label="",title=str)
#    StatsPlots.png(string(path,strname,"violin_all.png"))



    data11_main=Data[11][end-2:end]
    data15=Data[15][:,end-2:end]


    data11=log.(data11_main)./log(2)
    data11=repeat(data11',M,1)
    data11=reshape(data11,1,size(data11,1)*size(data11,2))

    data15=reshape(data15,size(data15,1)*size(data15,2),1)
    gr(size = (1000, 1000))
    df = DataFrame(sample = vec(data11), res= vec(data15))
    @df df violin(:sample,:res,linewidth=0,xticks=log.(Data[11])./log(2),xlabel="log2 of samples",ylabel="expected value",label="")
    @df df boxplot!(:sample,:res,linewidth=1.5,label="",fillalpha=0.75)
    @df df dotplot!(:sample,:res,markersize = 2.5,markerstrokewidth = 1.5,markerstrokealpha = 0.2, markerstrokecolor = :black,xticks=log.(data11_main)./log(2),label="",title=str)
    StatsPlots.png(string(path,strname,"violin_end.png"))






    data11=Data[11]
    data15=0
    GC.gc()
    data15=Data[15]

    sol = analytical_sol(1000,s,Data[13])
    println(sol)
    solMat=ones(length(data15),1)*sol


    data11=log.(Data[11])./log(2)
    data11=repeat(data11',M,1)
    data11=reshape(data11,1,size(data11,1)*size(data11,2))

    data15=(abs.(reshape(data15,size(data15,1)*size(data15,2),1).-solMat) ./solMat )
    data15[data15 .== 0] .=  minimum(data15[data15 .!= 0])
    data15=log.(data15)/log(10)


    mn = mean((Data[15] .- sol)./sol,dims=1)



    gr(size = (1000, 1000))
    df = DataFrame(sample = vec(data11), res= vec(data15))

    @df df violin(:sample,:res,linewidth=0,xticks=Data[11],xlabel="samples in base 2",ylabel="rel abs error in base 10",label="")
    @df df boxplot!(:sample,:res,linewidth=1.5,label="",fillalpha=0.75,yticks=collect(0:-1:-16))
    @df df dotplot!(:sample,:res,markersize = 2.5,markerstrokewidth = 1.5,markerstrokealpha = 0.2, markerstrokecolor = :black,label="",title=str,xticks=log.(Data[11])./log(2),yticks=collect(0:-1:-16),ylims=(-14,0))
#    StatsPlots.png(string(path,strname,"violin_error.png"))


end




main()
