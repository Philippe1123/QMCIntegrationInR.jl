using DigitalNets
using LatticeRules

#using PyPlot

using Statistics: mean, std
using SpecialFunctions: erf, erfinv, gamma, gamma_inc
using StringLiterals
using PrettyTables
using Random
using DelimitedFiles
using JLD2
using FileIO
using PyPlot
using FiniteElementDiffusion



NormalPdf(x) = 1 / sqrt(2 * pi) * exp(-1 / 2 * x^2)



function main()


    s = 3
    M = 32


    tol = 1 * 10.0 .^ -9 .* 3 .^ (20 .- 1 .- (1.0:1.0:19.0))
    #    tol = [3*10.0 .^ -8]

    Data = RunSimulation(
        s,
        M,
        tol,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    #writeOut(Data, "1_v8_res")

end




function RunSimulation(
    s::Int64,
    M::Int64,
    tol::Vector,
    QMCGenerator::Union{DigitalNet64,LatticeRule},
)
    alpha = 3 #hardcode alpha



    estimatedTime = []
    estimatedTruncationErrors = []
    estimatedCubatureErrors = []







    BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
    SampleExponent = 2


    for tolerance in tol
        SampleExponent = 2
        println("Currently running tolerance number ", tolerance)
        truncationError = 10 # initalisation
        QMCResultsInternals = [] # reintilized for each tolerance, only used internally
        estimatedCubatureErrorsInternals = [] # reintilized for each tolerance, only used internally
        estimatedTruncationErrorsInternals = [] # reintilized for each tolerance, only used internally
        samplesInternals = [] # reintilized for each tolerance, only used internally
        boundsOfBoxesInternals = [] # reintilized for each tolerance, only used internally
        t = @elapsed begin
            ########################################## Adjust truncation error error
            while truncationError > tolerance / 2
                BoxBoundary = 0
                numberOfPointsBox = 2^(SampleExponent)
                push!(samplesInternals, numberOfPointsBox)
                pointsBox = zeros(s, numberOfPointsBox, M)

                # We first generate *all* the points for all the shifts...
                BoxBoundary = sqrt(2 * alpha * log(numberOfPointsBox))

                push!(boundsOfBoxesInternals, BoxBoundary)
                Random.seed!(1234)
                for shiftId = 1:M
                    shiftedQMCGenerator = randomizedGenerator(QMCGenerator)
                    for id = 1:numberOfPointsBox
                        pointsBox[:, id, shiftId] =
                            map.(
                                x -> -BoxBoundary + (BoxBoundary - (-BoxBoundary)) * x,
                                shiftedQMCGenerator[id-1],
                            )
                    end
                end

                # Solve the problem
                G_fine = SolveRoutine(pointsBox)


                QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                push!(QMCResultsInternals, QMC_Q[1])
                push!(estimatedCubatureErrorsInternals, QMC_std)
                pointsBox = 0
                G_fine = 0
                GC.gc()


                # Routine to check evolution of truncation error
                if length(QMCResultsInternals) > 1
                    println(QMCResultsInternals[end-1])
                    println(QMCResultsInternals[end])
                    truncationError =
                        abs.(QMCResultsInternals[end-1] - QMCResultsInternals[end]) /
                        QMCResultsInternals[end]
                    push!(estimatedTruncationErrorsInternals, truncationError)
                else
                    push!(estimatedTruncationErrorsInternals, -1)
                end





                SampleExponent = SampleExponent + 1
            end




            ####################### Adjust Cubature error
            estimatedCubatureError = estimatedCubatureErrorsInternals[end]
            SampleExponent = 2
            estimatedCubatureError = 10

            while estimatedCubatureError > tolerance / 2
                BoxBoundary = boundsOfBoxesInternals[end]
                numberOfPointsBox = 2^(SampleExponent)
                pointsBox = zeros(s, numberOfPointsBox, M)
                Random.seed!(1234)
                for shiftId = 1:M
                    shiftedQMCGenerator = randomizedGenerator(QMCGenerator)
                    for id = 1:numberOfPointsBox
                        pointsBox[:, id, shiftId] =
                            map.(
                                x -> -BoxBoundary + (BoxBoundary - (-BoxBoundary)) * x,
                                shiftedQMCGenerator[id-1],
                            )
                    end
                end

                G_fine = SolveRoutine(pointsBox)



                # Compute std of qmc and expected value
                QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                estimatedCubatureError = QMC_std
                push!(estimatedCubatureErrorsInternals, QMC_std)
                pointsBox = 0
                G_fine = 0
                GC.gc()

                SampleExponent = SampleExponent + 1

            end
 
        end # end of @elapsed

        push!(estimatedCubatureErrors, estimatedCubatureErrorsInternals[end])
        push!(estimatedTruncationErrors, estimatedTruncationErrorsInternals[end])
        push!(estimatedTime, t)

        println(
            t,
            "  ",
            estimatedCubatureErrorsInternals[end],
            "  ",
            estimatedTruncationErrorsInternals[end],
        )



    end

    figure()
    loglog(estimatedTime, estimatedTruncationErrors, "*-k")
    loglog(estimatedTime, estimatedCubatureErrors, "*-r")

    loglog(estimatedTime, estimatedTime .^ -1, "--b")
    loglog(estimatedTime, estimatedTime .^ -2, "--r")
    loglog(estimatedTime, estimatedTime .^ -3, "--y")
    grid(which = "both", ls = "-")

    """
    println(QMCType)
    println("alpha is ", alpha)
    println("exact solution is ", exactSol)
    println("stochastic dimensions is equal to ", s)
    println("params is equal to ", params)
    println("number of shifts is ", M)


    data = hcat(
        log.(N) ./ log(2),
        N,
        boundsOfBoxes,
        QMCResults,
        trueErrors,
        solutionsOnBox,
        trueTruncationErrors,
        trueCubatureErrors,
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
    ])

    formatters = ft_printf(
        ["%-3d", "%-3d", "%16.8f", "%16.16f", "%.5e", "%16.16f", "%.5e", "%.5e"],
        [1, 2, 3, 4, 5, 6, 7, 8],
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
    """
end


function SolveRoutine(pointsBox::Array)
    G_fine = prod(
        (1 .+ abs.(pointsBox) .^ 2.6) .* 1 / (sqrt(2 * pi)) .* exp.(-(pointsBox .^ 2) ./ 2),
        dims = 1,
    )
    return G_fine

end

function ComputeQMCStdAndExp(G_fine::Array, BoxBoundary::Float64, s::Int64, M::Int64)


    # pointsBox is s-by-N-by-M --f--> 1-by-N-by-M
    # G_fine = mean(f(pointsBox, params), dims = 2) # 1-by-1-by-M
    G_fine = mean(G_fine, dims = 2)
    QMC_R = abs.(G_fine) * (BoxBoundary * 2)^s
    QMC_Q = mean(QMC_R, dims = 3)
    QMC_std = std(QMC_R) / sqrt(M)
    return QMC_std, QMC_Q

end





main()
