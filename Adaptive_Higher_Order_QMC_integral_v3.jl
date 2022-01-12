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



NormalPdf(x) = 1 / sqrt(2 * pi) * exp(-1 / 2 * x^2)


analytical_sol(a::Real, s::Int64, sigma::Real) =
    (
        (
            gamma((1 + sigma) / 2) -
            gamma((1 + sigma) / 2) * gamma_inc((1 + sigma) / 2, a^2 / 2, 0)[2]
        ) * 2^(sigma / 2) / sqrt(pi) + erf(a / sqrt(2))
    )^s



function main()


    s = 3
    M = 32


    tol = 10.0 .^ (-1:-1:-9)

    Data = RunSimulation(
        s,
        M,
        tol,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )

end

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
    DictOfEstimatedTruncationErrors = Dict()
    DictOfEstimatedCubatureErrors = Dict()
    DictOfEstimatedCubatureErrorsTimings = Dict()


    exactCubatureErrors = []
    exactTruncationErrors = []





    BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
    SampleExponentBox = 2
    counter = 1
    corrfactor = 1

    for tolerance in tol
        #SampleExponentBox = 2
        println("##############################################")
        println("Currently running tolerance number ", tolerance)
        println("")
        truncationError = 10 * tolerance # initalisation
        QMCResultsInternals = [] # reintilized for each tolerance, only used internally
        estimatedCubatureErrorsInternals = [] # reintilized for each tolerance, only used internally
        estimatedCubatureErrorsInternalsTimings = []
        estimatedTruncationErrorsInternals = [] # reintilized for each tolerance, only used internally
        samplesInternals = [] # reintilized for each tolerance, only used internally
        boundsOfBoxesInternals = [] # reintilized for each tolerance, only used internally

        if length(estimatedTruncationErrors) > 0 &&
           estimatedTruncationErrors[end] < tolerance / 2 &&
           estimatedCubatureErrors[end] < tolerance / 2


            estimatedTruncationErrorsInternals = DictOfEstimatedTruncationErrors[counter-1]
            estimatedCubatureErrorsInternals = DictOfEstimatedCubatureErrors[counter-1]
            estimatedCubatureErrorsInternalsTimings =
                DictOfEstimatedCubatureErrorsTimings[counter-1]

            t = estimatedTime[end]
            println("Truncation and Cubature were already satisfied in previous run")


        else



            t = @elapsed begin
                ########################################## Adjust truncation error error
                while truncationError > tolerance / 2

                    timingQMC = @elapsed begin
                        BoxBoundary = 0
                        numberOfPointsBox = 2^(SampleExponentBox)
                        push!(samplesInternals, numberOfPointsBox)

                        #compute the box boundary
                        BoxBoundary = sqrt(2 * alpha * log(numberOfPointsBox))
                        push!(boundsOfBoxesInternals, BoxBoundary)
                        numberOfPointsBox = Int64(numberOfPointsBox * corrfactor) # apply correction after computation box

                        # compute all the points
                        pointsBox = mapPoints(
                            M,
                            QMCGenerator,
                            numberOfPointsBox,
                            s,
                            BoxBoundary,
                        )

                        # Solve the problem
                        G_fine = SolveRoutine(pointsBox)
                        QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                        QMC_Q_final = QMC_Q[1]

                    end
                    push!(estimatedCubatureErrorsInternals, QMC_std)

                    if length(estimatedCubatureErrorsInternalsTimings) > 0
                        push!(
                            estimatedCubatureErrorsInternalsTimings,
                            timingQMC + estimatedCubatureErrorsInternalsTimings[end],
                        )
                    else
                        push!(estimatedCubatureErrorsInternalsTimings, timingQMC)
                    end




                    ####################### Adjust Cubature error
                    estimatedCubatureError = QMC_std
                    SampleExponentCubature = SampleExponentBox    # start with exponent used for box
                    searched = false

                    #########################################################################################
                    while estimatedCubatureError > tolerance / 2 
                        timingQMC = @elapsed begin

                            numberOfPointsBox = Int64(2^(SampleExponentCubature) * corrfactor)
                            pointsBox = mapPoints(
                                M,
                                QMCGenerator,
                                numberOfPointsBox,
                                s,
                                BoxBoundary,
                            )
                            G_fine = SolveRoutine(pointsBox)
                            # Compute std of qmc and expected value
                            QMC_std, QMC_Q =
                                ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                            estimatedCubatureError = QMC_std
                            push!(estimatedCubatureErrorsInternals, QMC_std)

                        end
                        QMC_Q_final = QMC_Q[1]

                        #println(length(estimatedCubatureErrorsInternalsTimings))
                        push!(
                            estimatedCubatureErrorsInternalsTimings,
                            timingQMC + estimatedCubatureErrorsInternalsTimings[end],
                        )
                        println(
                            "qmc error is ",
                            QMC_std,
                            " on box ",
                            -BoxBoundary,
                            " ",
                            BoxBoundary,
                            " exp val is ",
                            QMC_Q_final,
                            " eact sol is ",
                            analytical_sol(BoxBoundary, s, 2.6),
                        )
                        if estimatedCubatureError > tolerance / 2
                            corrfactor = corrfactor * 2
                        end
                        searched = true

                    end


                    while estimatedCubatureError < tolerance / 2 && searched == false
                        timingQMC = @elapsed begin

                            numberOfPointsBox = Int64(2^(SampleExponentCubature) * corrfactor)
                            println(numberOfPointsBox)
                            pointsBox = mapPoints(
                                M,
                                QMCGenerator,
                                numberOfPointsBox,
                                s,
                                BoxBoundary,
                            )
                            G_fine = SolveRoutine(pointsBox)
                            # Compute std of qmc and expected value
                            QMC_std, QMC_Q =
                                ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                            estimatedCubatureError = QMC_std

                        end
                        #println(length(estimatedCubatureErrorsInternalsTimings))

                        if estimatedCubatureError < tolerance / 2
                            searched = true 
                            QMC_Q_final = QMC_Q[1]
                            corrfactor = corrfactor / 2
                            push!(
                                estimatedCubatureErrorsInternalsTimings,
                                timingQMC + estimatedCubatureErrorsInternalsTimings[end],
                            )
                            push!(estimatedCubatureErrorsInternals, QMC_std)

                            println(
                                "qmc error is ",
                                QMC_std,
                                " on box ",
                                -BoxBoundary,
                                " ",
                                BoxBoundary,
                                " exp val is ",
                                QMC_Q_final,
                                " eact sol is ",
                                analytical_sol(BoxBoundary, s, 2.6),
                            )
                        end
                    end
                    #########################################################################################





                    println(
                            "qmc error is ",
                            QMC_std,
                            " on box ",
                            -BoxBoundary,
                            " ",
                            BoxBoundary,
                            " exp val is ",
                            QMC_Q_final,
                            " eact sol is ",
                            analytical_sol(BoxBoundary, s, 2.6),
                        )
                    println("------Finished Cubature error--------")

                    push!(QMCResultsInternals, QMC_Q_final)

                    # Routine to check evolution of truncation error
                    if length(QMCResultsInternals) > 1
                        truncationError =
                            abs.(QMCResultsInternals[end-1] - QMCResultsInternals[end]) / QMCResultsInternals[end]
                        println("truncation error is ", truncationError)

                        push!(estimatedTruncationErrorsInternals, truncationError)
                    else
                        push!(estimatedTruncationErrorsInternals, -1)
                    end

                    if truncationError > tolerance / 2
                        estimatedCubatureErrorsInternals = [] # clear array when computing new truncation error
                        lastTiming = estimatedCubatureErrorsInternalsTimings[end]
                        estimatedCubatureErrorsInternalsTimings = []
                        push!(estimatedCubatureErrorsInternalsTimings, lastTiming)

                    end


                    SampleExponentBox = SampleExponentBox + 1 # increase exponent for next box

                end
                SampleExponentBox = SampleExponentBox - 1   #if while loop of truncation error is satisfied reset the sample exponent back with one, avoid unecessary increasing the starting box for new tolerances 


                println("##############################################")


            end # end of @elapsed
        end # end of if else



        #### computation of exact sol and updating arrays for printing

        DictOfEstimatedTruncationErrors[counter] = estimatedTruncationErrorsInternals
        DictOfEstimatedCubatureErrors[counter] = estimatedCubatureErrorsInternals
        DictOfEstimatedCubatureErrorsTimings[counter] =
            estimatedCubatureErrorsInternalsTimings



        exactSolOnBox = analytical_sol(BoxBoundary, s, 2.6)

        exactSol = analytical_sol(1000, s, 2.6)

        push!(
            exactCubatureErrors,
            abs(exactSolOnBox - QMCResultsInternals[end]) ./ QMCResultsInternals[end],
        )
        push!(exactTruncationErrors, abs(exactSol - exactSolOnBox) ./ exactSol)

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


        counter = counter + 1
    end

    figure()
    loglog(estimatedTime, estimatedTruncationErrors, "*-k")
    loglog(estimatedTime, estimatedCubatureErrors, "*-r")





    loglog(estimatedTime, estimatedTime .^ -1, "--b")
    loglog(estimatedTime, estimatedTime .^ -2, "--r")
    loglog(estimatedTime, estimatedTime .^ -3, "--y")
    loglog(estimatedTime, exactTruncationErrors, "*--k")
    loglog(estimatedTime, exactCubatureErrors, "*--r")
    grid(which = "both", ls = "-")
    legend((
        "estimatedTruncationErrors",
        "estimatedCubatureErrors",
        "time^-1",
        "time^-2",
        "time^-3",
        "exact truncation error",
        "exact cubature error",
    ))
    xlabel("time [sec]")
    ylabel("error [/]")
    for i = 1:length(tol)


        loglog(
            DictOfEstimatedCubatureErrorsTimings[i][2:end],
            DictOfEstimatedCubatureErrors[i][1:end],
            "-.or",
            alpha = 0.3,
            mec = "r",
        )
    end
    #    println(DictOfEstimatedCubatureErrorsTimings[1])
    #    println(DictOfEstimatedCubatureErrorsTimings[1][end])
    #    println(DictOfEstimatedCubatureErrorsTimings[2][end])
    #    println(estimatedTime)

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


function mapPoints(
    M::Int64,
    QMCGenerator::Union{DigitalNet64,LatticeRule},
    numberOfPointsBox::Int64,
    s::Int64,
    BoxBoundary::Float64,
)

    Random.seed!(1234)
    pointsBox = zeros(s, numberOfPointsBox, M)
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
    return pointsBox
end






main()
