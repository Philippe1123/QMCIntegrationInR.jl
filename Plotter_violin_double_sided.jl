
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

    for id = 13:2:17
        Data_1 = load(string(@__DIR__,"/res/NoMapping/31082021/",id, ".jld2"))["data"]
        Data_2 = load(string(@__DIR__,"/res/NoMapping/31082021/",id+1, ".jld2"))["data"]
        plotter(Data_1,Data_2)

    end


end


function plotter(Data_1::Dict,Data_2::Dict)


    QMCType1 = Data_1[8]
    s1 = Data_1[9]
    M1 = Data_1[10]
    nbOfSamples1 = Data_1[11]



    QMCType2 = Data_2[8]
    s2 = Data_2[9]
    M2 = Data_2[10]
    nbOfSamples2 = Data_2[11]


   


    str1 = string(
        QMCType1,
        "\n",
        "s = ",
        s1,
        ", shift = ",
        M1,
        ", params = ",
        Data_1[13],
        ", alpha = ",
        Data_1[14],
    )


    str2 = string(
        QMCType2,
        "\n",
        "s = ",
        s2,
        ", shift = ",
        M2,
        ", params = ",
        Data_2[13],
        ", alpha = ",
        Data_2[14],
    )


    strname=string(
        QMCType1,"_",QMCType2,
        "_",
        "s=",
        s1,
        "_shift=",
        M1,
        "_params=",
        Data_1[13],
        "_alpha=",
        Data_1[14],
    )

"""
    data11=Data_1[11]
    data15=Data_1[15]
    data11_plot=Data_1[11]

    data11=log.(data11)./log(2)
    data11=repeat(data11',M1,1)
    data11=reshape(data11,1,size(data11,1)*size(data11,2))
    data15=reshape(data15,size(data15,1)*size(data15,2),1)
    gr(size = (1000, 1000))
    df = DataFrame(sample = vec(data11), res= vec(data15))
    @df df violin(:sample,:res,linewidth=0,xticks=log.(data11_plot)./log(2),xlabel="log2 of samples",ylabel="expected value",label=str1,side=:left)
    #@df df boxplot!(:sample,:res,linewidth=1.5,label="",fillalpha=0.75,side=:left)
    @df df dotplot!(:sample,:res,markersize = 2.5,markerstrokewidth = 1.5,markerstrokealpha = 0.2, markerstrokecolor = :black,xticks=log.(data11_plot)./log(2),label="",side=:left)
    StatsPlots.png(string(strname,"violin_all.png"))



    data11_main=Data_1[11][end-2:end]
    data15=Data_1[15][:,end-2:end]
    data11=log.(data11_main)./log(2)
    data11=repeat(data11',M1,1)
    data11=reshape(data11,1,size(data11,1)*size(data11,2))

    data15=reshape(data15,size(data15,1)*size(data15,2),1)
    gr(size = (1000, 1000))
    df = DataFrame(sample = vec(data11), res= vec(data15))
    @df df violin(:sample,:res,linewidth=0,xticks=log.(data11_main)./log(2),xlabel="log2 of samples",ylabel="expected value",label=str1,side=:left)
    #@df df boxplot!(:sample,:res,linewidth=1.5,label="",fillalpha=0.75,side=:left)
    @df df dotplot!(:sample,:res,markersize = 2.5,markerstrokewidth = 1.5,markerstrokealpha = 0.2, markerstrokecolor = :black,xticks=log.(data11_main)./log(2),label="",side=:left)
    StatsPlots.png(string(strname,"violin_end.png"))

"""


 
#############

    lnData1=length(Data_1[11])
    lnData2=length(Data_2[11])



    data11=Data_1[11]
    data15=Data_1[15]
    data11_plot=Data_1[11]


    sol = analytical_sol(1000,s1,Data_1[13])
    solMat=ones(length(data15),1)*sol
    data11=log.(Data_1[11])./log(2)
    data11=repeat(data11',M1,1)
    data11=reshape(data11,1,size(data11,1)*size(data11,2))
    data15=(  log.(abs.(reshape(data15,size(data15,1)*size(data15,2),1).-solMat) ./solMat ))/log(10)
    gr(size = (1000, 1000))
    df = DataFrame(sample = vec(data11), res= vec(data15))
    @df df violin(:sample,:res,linewidth=0,xticks=data11_plot,xlabel="samples in base 2",ylabel="rel abs error in base 10",label=str1,side=:left)
    #@df df boxplot!(:sample,:res,linewidth=1.5,label="",fillalpha=0.75,yticks=collect(0:-1:-16),side=:left)
    @df df dotplot!(:sample,:res,markersize = 1.5,markerstrokewidth = 1.5,markerstrokealpha = 0.2, markerstrokecolor = :black,xticks=log.(data11_plot)./log(2),yticks=collect(0:-1:-16),ylims=(-14,0),label="",side=:left)
    
    data11=Data_2[11][1:lnData1]
    data15=Data_2[15][:,1:lnData1]
    data11_plot=Data_1[11][1:lnData1]


    sol = analytical_sol(1000,s2,Data_2[13])
    solMat=ones(length(data15),1)*sol
    data11=log.(Data_1[11])./log(2)
    data11=repeat(data11',M1,1)
    data11=reshape(data11,1,size(data11,1)*size(data11,2))
    data15=(  log.(abs.(reshape(data15,size(data15,1)*size(data15,2),1).-solMat) ./solMat ))/log(10)
    gr(size = (1000, 1000))
    df = DataFrame(sample = vec(data11), res= vec(data15))
    @df df violin!(:sample,:res,linewidth=0,xticks=data11_plot,xlabel="samples in base 2",ylabel="rel abs error in base 10",label=str2,side=:right)
    #@df df boxplot!(:sample,:res,linewidth=1.5,label="",fillalpha=0.75,yticks=collect(0:-1:-16),side=:left)
    @df df dotplot!(:sample,:res,markersize = 1.5,markerstrokewidth = 1.5,markerstrokealpha = 0.2, markerstrokecolor = :black,xticks=log.(data11_plot)./log(2),yticks=collect(0:-1:-16),ylims=(-14,0),label="",side=:right)
    
    
 #############  
    
    StatsPlots.png(string(strname,"violin_error.png"))





end




main()
