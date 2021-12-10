using DigitalNets
using LatticeRules


using Statistics: mean, std
using SpecialFunctions: erf, erfinv, gamma, gamma_inc
using StringLiterals
using PrettyTables
using Random
using DelimitedFiles
using JLD2
using FileIO
using PyPlot

Φ⁻¹(x::T where {T<:Real}) = √2 * erfinv(2 * x - 1)
Φ⁻¹(x::T where {T<:Real},σ::T where {T<:Real}) = √2 * erfinv(2 * x - 1) * σ

Φ(x::T where {T<:Real}) = 1 / 2 * (1 + erf(x / √2))
Φ(x::T where {T<:Real},σ::T where {T<:Real}) = 1 / 2 * (1 + erf(x /(σ * √2)))
cdfNorm(x::T where {T<:Real}) = 1 / (sqrt(2 * pi)) .* exp.(-(x .^ 2) ./ 2)



#analytical_sol(a::Real,s::Int64) = ((gamma(4/5)/(2^(1/5))-gamma(4/5)*gamma_inc(4/5,a^2/2,0)[2]/(2^(1/5)))*2/sqrt(2*pi)+(Φ(a) - Φ(-a)))^s
analytical_sol(a::Real, s::Int64, sigma::Real) =
    (
        (
            gamma((1 + sigma) / 2) -
            gamma((1 + sigma) / 2) * gamma_inc((1 + sigma) / 2, a^2 / 2, 0)[2]
        ) * 2^(sigma / 2) / sqrt(pi) + erf(a / sqrt(2))
    )^s



    function main()

        #### Input parameters
        M = 64 # number of shifts
        N_lattice = 2 .^ collect(4:1:20)
        N_net = 2 .^ collect(4:1:20)


        #  generator = DigitalNet64InterlacedTwo(s)
        #generator = DigitalNet64InterlacedThree(s)

        #generator = LatticeRule(s)

        #generator = DigitalNet64(s)


        # dim = 1

        s = 1


        Data = RunSimulation(
            s,
            M,
            N_lattice,
            1,
            0.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "1_v6_res")

        Data = RunSimulation(s, M, N_net, 1, 0.6, DigitalNet64(s),true)
        writeOut(Data, "2_v6_res")

        Data = RunSimulation(
            s,
            M,
            N_lattice,
            2,
            1.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "3_v6_res")

        Data = RunSimulation(s, M, N_net, 2, 1.6, DigitalNet64InterlacedTwo(s),true)
        writeOut(Data, "4_v6_res")
        Data = RunSimulation(
            s,
            M,
            N_lattice,
            3,
            2.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "5_v6_res")

        Data = RunSimulation(s, M, N_net, 3, 2.6, DigitalNet64InterlacedThree(s),true)
        writeOut(Data, "6_v6_res")

"""
        # dim = 2
        s = 2
        Data = RunSimulation(
            s,
            M,
            N_lattice,
            1,
            0.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "7_Inv")

        Data = RunSimulation(s, M, N_net, 1, 0.6, DigitalNet64(s),true)
        writeOut(Data, "8_Inv")


        Data = RunSimulation(
            s,
            M,
            N_lattice,
            2,
            1.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "9_Inv")

        Data = RunSimulation(s, M, N_net, 2, 1.6, DigitalNet64InterlacedTwo(s),true)
        writeOut(Data, "10_Inv")

        Data = RunSimulation(
            s,
            M,
            N_lattice,
            3,
            2.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "11_Inv")

        Data = RunSimulation(s, M, N_net, 3, 2.6, DigitalNet64InterlacedThree(s),true)
        writeOut(Data, "12_Inv")

        # dim = 3
        s = 3

        Data = RunSimulation(
            s,
            M,
            N_lattice,
            1,
            0.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "13_Inv")
        Data = RunSimulation(s, M, N_net, 1, 0.6, DigitalNet64(s),true)
        writeOut(Data, "14_Inv")
        Data = RunSimulation(
            s,
            M,
            N_lattice,
            2,
            1.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "15_Inv")

        Data = RunSimulation(s, M, N_net, 2, 1.6, DigitalNet64InterlacedTwo(s),true)
        writeOut(Data, "16_Inv")

        Data = RunSimulation(
            s,
            M,
            N_lattice,
            3,
            2.6,
            LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),true
        )
        writeOut(Data, "17_inv")

        Data = RunSimulation(s, M, N_net, 3, 2.6, DigitalNet64InterlacedThree(s),true)
        writeOut(Data, "18_Inv")
"""


    end # end of function main()


## Quick and dirty solution for the "isa" problem (needs to be fixed in a more decent matter)
"""
Create a randomized generator with a random shift or random digital shift for the passed in QMC generator.
    """
randomizedGenerator(latticeGenerator::LatticeRule) = ShiftedLatticeRule(latticeGenerator)
function randomizedGenerator(digitalnetGenerator::DigitalNet64)
    DigitalShiftedDigitalNets64(digitalnetGenerator)
end
"""
String label for the current generator, either Lattice or Sobol1, Sobol2, or Sobol3.
    """
labelThisGenerator(latticeGenerator::LatticeRule) = "Lattice"
labelThisGenerator(digitalnetGenerator::DigitalNet64) =
    "Sobol$(Int32(round(sqrt(reversebits(digitalnetGenerator.C[1])))))"


"""
RunSimulation(s, M, requestedTolerances, QMCGenerator)

Run the algorithm for our test function in `s` dimensions with `M` shifts and a
    starting number of points `2^N0` for all the requested tolerances in
    `requestedTolerances` using the QMC point generator in `QMCGenerator`.
    """
function RunSimulation(
    s::Int64,
    M::Int64,
    N::Vector,
    alpha::Int64,
    params::Float64,
    QMCGenerator::Union{DigitalNet64,LatticeRule},
    correctionFactor::Bool,
)

    QMCType = labelThisGenerator(QMCGenerator)



    timings = zeros(length(N))
    trueErrors = zeros(length(N))
    estimatedTruncationErrors = zeros(length(N))
    estimatedCubatureErrors = zeros(length(N))
    trueTruncationErrors = zeros(length(N))
    trueCubatureErrors = zeros(length(N))
    nbOfSamples = zeros(length(N))
    maxLevelPerRequestedTol = Int64.(zeros(length(N)))
    boundsOfBoxes = zeros(length(N))
    QMCResults = zeros(length(N))
    solutionsOnBox = zeros(length(N))
    correctionFactors = zeros(length(N))
    exactSol = 0
    shiftAverages = zeros(M, length(N))

    f(x, params, a) = prod((1 .+ abs.(a*x) .^ params) .* a .* exp.(-x.^2 ./ 2 .* (a^2-1)), dims = 1)


    # why do these two variables need to be declared here? can't we rewrite the logic?
    BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
    #truncationError = 10000
    #cubatureError = 10000
    soll = 10000
    idx = 1
    for ell in N
        println("Currently running sample number ", ell, " exponent ", log.(ell) ./ log(2))

        t = @elapsed begin

            cubature_error = 1000 # initalisation
            truncation_error = 1000 # initalisation
            QMC_Q = 0
            corrfactor_fine = 0
            BoxBoundary = 0
            exactTruncationError = 1000 # initalisation
            exactCubatureError = 1000 # initalisation
            totError = 1000


            numberOfPointsBox = ell
            pointsBox = zeros(s, numberOfPointsBox, M)
            pointsBox_unmapped = zeros(s, numberOfPointsBox, M)

            # We first generate *all* the points for all the shifts...
            BoxBoundary = sqrt(2 * alpha * log(numberOfPointsBox))
            boundsOfBoxes[idx] = BoxBoundary

            Random.seed!(1234)
            for shiftId = 1:M
                shiftedQMCGenerator = randomizedGenerator(QMCGenerator)

                for id = 1:numberOfPointsBox
                    pointsBox[:, id, shiftId] =
                        map.(
                            x -> Φ⁻¹(
                                Φ(-BoxBoundary,1) + (Φ(BoxBoundary,1) - Φ(-BoxBoundary,1)) * x,1.
                            ),
                            shiftedQMCGenerator[id-1],
                        )
                        pointsBox_unmapped[:, id, shiftId]=shiftedQMCGenerator[id-1]


                end
            end




            # pointsBox is s-by-N-by-M --f--> 1-by-N-by-M
            G_fine = mean(f(pointsBox,params,10), dims = 2) # 1-by-1-by-M
#            pointsBox = 0
            GC.gc()

 #           println(pointsBox)
 #           sleep(1.0)
 """
            figure()
            res=f(pointsBox,params)
            for ip=1:size(pointsBox,2)
                plot(pointsBox[1,ip, 1],res[1,ip,1],"*k")
            end

            figure()
            for ip=1:size(pointsBox,2)
                plot(pointsBox[1,ip, 1],"*k")
                sleep(0.01)
            end
            """

            corrfactor_fine =
                correctionFactor == true ? (Φ(BoxBoundary,1.) - Φ(-BoxBoundary,1.))^s : 1

            QMC_R = abs.(G_fine) * (corrfactor_fine)
            G_fine = 0
            GC.gc()
            shiftAverages[:, idx] = vec(reshape(QMC_R, M, 1, 1))


            QMC_Q = mean(QMC_R, dims = 3)
            QMCResults[idx] = QMC_Q[1]


            QMC_std = std(QMC_R) / sqrt(M)
            QMC_R = 0
            GC.gc()



            cubature_error = QMC_std





            exactSolOnBox = analytical_sol(BoxBoundary, s, params)
            solutionsOnBox[idx] = exactSolOnBox
            exactSol = analytical_sol(1000, s, params)
            exactCubatureError = abs(exactSolOnBox - QMC_Q[1]) ./ QMC_Q[1]
            exactTruncationError = abs(exactSol - exactSolOnBox) ./ exactSol
            totError = abs(exactSol - QMC_Q[1]) ./ exactSol






        end # end of @elapsed


        correctionFactors[idx] = corrfactor_fine
        estimatedTruncationErrors[idx] = truncation_error
        estimatedCubatureErrors[idx] = cubature_error
        timings[idx] = t
        trueErrors[idx] = totError
        trueTruncationErrors[idx] = exactTruncationError
        trueCubatureErrors[idx] = exactCubatureError
        idx = idx + 1

    end

    println(QMCType)
    println("alpha is ", alpha)
    println("exact solution is ", exactSol)
    println("stochastic dimensions is equal to ", s)
    println("params is equal to ", params)
    println("number of shifts is ", M)
    println("Correction factor enabled ?", correctionFactor)


    data = hcat(
        log.(N) ./ log(2),
        N,
        boundsOfBoxes,
        QMCResults,
        trueErrors,
        solutionsOnBox,
        trueTruncationErrors,
        trueCubatureErrors,
        correctionFactors,
    )
    header = ([
        "m",
        "n",
        "a",
        "Q",
        "tot rel error",
        "Iab",
        "trunc rel error (exact)",
        "box cub rel error (exact)",
        "Corr. fact"
    ])

    formatters = ft_printf(
        ["%-3d", "%-3d", "%16.8f", "%16.16f", "%.5e", "%16.16f", "%.5e", "%.5e","%16.8f"],
        [1, 2, 3, 4, 5, 6, 7, 8,9],
    )

    pretty_table(data; header = header, formatters = formatters)




    Data = Dict()
    Data[2] = timings
    Data[3] = trueErrors
    Data[4] = estimatedCubatureErrors
    Data[5] = trueCubatureErrors
    Data[6] = estimatedTruncationErrors
    Data[7] = trueTruncationErrors
    Data[8] = QMCType
    Data[9] = s
    Data[10] = M
    Data[11] = N
    Data[13] = params
    Data[14] = alpha
    Data[15] = shiftAverages

    return Data
end



function plotter(Data::Dict)

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
    println(timings)
    #loglog(nbOfSamples, abs.(estimatedCubatureErrors), "-r*")
    loglog(nbOfSamples, abs.(trueCubatureErrors), "-rs")
    #loglog(nbOfSamples, abs.(estimatedTruncationErrors), "-g*")
    loglog(nbOfSamples, abs.(trueTruncationErrors), "-gs")
    loglog(nbOfSamples, nbOfSamples .^ -1, "--b")
    loglog(nbOfSamples, nbOfSamples .^ -2, "--r")
    loglog(nbOfSamples, nbOfSamples .^ -3, "--y")


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
        Data[14]
    )
    title(str)
    ylabel("Abs. error")
    xlabel("N")
    legend((
        "true error",
        "true cubature error",
        "true truncation error",
        "N^-1",
        "N^-2",
        "N^-3",
    ))
    println(estimatedTruncationErrors)
    println(estimatedCubatureErrors)
    println(trueTruncationErrors)
    println(trueCubatureErrors)
    savefig(string(str,".png"))


end


function writeOut(Data::Dict, str::String)

    currentFolder = @__DIR__
    Path = string(currentFolder, "/", str, ".jld2")
    println(Path)
    save(Path, "data", Data)


end

main()
