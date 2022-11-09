
###########THIS IS THE LATEST VERSION FOR THE INTEGRAL PROBLEM
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
using Distributed
using LaTeXStrings

using DigitalNets
using LatticeRules




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


    tol = 10.0 .^ (-1:-1:-10)

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
    DictOfSamples = Dict()
    SampleExponentCubatureArray = []

    exactCubatureErrors = []
    exactTruncationErrors = []





    BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
    SampleExponentBox = 2
    SampleExponentCubature = SampleExponentBox
    counter = 1
    for tolerance in tol
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
        numberofSamples = []
        """    
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
      """


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
                    if (length(SampleExponentCubatureArray) == 0)
                        SampleExponentCubature = SampleExponentBox    # start with exponent used for box
                    else
                        SampleExponentCubature = SampleExponentCubatureArray[end]
                    end

                    numberOfPointsBox = 2^(SampleExponentCubature)

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
                end

                println(
                    "qmc error is ",
                    QMC_std,
                    " on box ",
                    -BoxBoundary,
                    " ",
                    BoxBoundary,
                    " exp val is ",
                    QMC_Q,
                    " exact sol is ",
                    analytical_sol(BoxBoundary, s, 2.6),
                    " samples is ",
                    numberOfPointsBox,
                )

                push!(estimatedCubatureErrorsInternals, QMC_std)

                if length(estimatedCubatureErrorsInternalsTimings) > 0
                    push!(
                        estimatedCubatureErrorsInternalsTimings,
                        timingQMC + estimatedCubatureErrorsInternalsTimings[end],
                    )
                else
                    push!(estimatedCubatureErrorsInternalsTimings, timingQMC)
                end



                push!(numberofSamples, numberOfPointsBox)

                ####################### Adjust Cubature error
                estimatedCubatureError = QMC_std
  
                timingsearch = 0
                f = (x, y) -> x + y

                timingsearch = @elapsed begin
                    println(SampleExponentCubature)
                    searchDirct, out = searchDirection(
                        M,
                        s,
                        BoxBoundary,
                        SampleExponentCubature,
                        tolerance / 2,
                        QMCGenerator,
                    )

                    isSearched = false
                    f = (x, y) -> -10
                    if searchDirct == 1
                        isSearched = true
                        estimatedCubatureError = out[2][1]
                        QMC_Q = out[2][2] ################################################
                        f = (x, y) -> Int64(x + 0)
                        cond =
                            (estimatedCubatureError, tolerance, isSearched, QMC_std_next) ->
                                estimatedCubatureError > tolerance / 2 ||
                                    isSearched == false
                    elseif searchDirct == 2
                        estimatedCubatureError = out[1][1]
                        QMC_Q = out[1][2]##################################################
                        SampleExponentCubature = SampleExponentCubature - 1
                        f = (x, y) -> max(Int64(x - y), 1)
                        cond =
                            (estimatedCubatureError, tolerance, isSearched, QMC_std_next) ->
                                (
                                    estimatedCubatureError < tolerance / 2 &&
                                    QMC_std_next < tolerance / 2 &&
                                    estimatedCubatureError == QMC_std_next
                                ) || isSearched == false
                    elseif searchDirct == 3
                        f = (x, y) -> Int64(x + y)
                        estimatedCubatureError = out[3][1]
                        QMC_Q = out[3][2]################################################################
                        SampleExponentCubature = SampleExponentCubature + 1
                        cond =
                            (estimatedCubatureError, tolerance, isSearched, QMC_std_next) ->
                                estimatedCubatureError > tolerance / 2 ||
                                    isSearched == false
                    end
                end
                estimatedCubatureErrorsInternalsTimings[end] =
                    estimatedCubatureErrorsInternalsTimings[end] + timingsearch

                QMC_std_next = estimatedCubatureError
                #### must force algorithm to search many times in other diretion, currently goes only once
                while cond(estimatedCubatureError, tolerance, isSearched, QMC_std_next)
                    isSearched = true
                    timingQMC = @elapsed begin
                        numberOfPointsBox = 2^(SampleExponentCubature)
                        println(numberOfPointsBox)

                        push!(numberofSamples, numberOfPointsBox)

                        pointsBox = mapPoints(
                            M,
                            QMCGenerator,
                            numberOfPointsBox,
                            s,
                            BoxBoundary,
                        )
                        G_fine = SolveRoutine(pointsBox)
                        # Compute std of qmc and expected value
                        QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                        estimatedCubatureError = QMC_std
                        push!(estimatedCubatureErrorsInternals, QMC_std)

                    end
                    push!(
                        estimatedCubatureErrorsInternalsTimings,
                        timingQMC + estimatedCubatureErrorsInternalsTimings[end],
                    )

                    if searchDirct == 2
                        timingQMC = @elapsed begin

                            numberOfPointsBox = 2^(SampleExponentCubature - 1)
                            pointsBox = mapPoints(
                                M,
                                QMCGenerator,
                                numberOfPointsBox,
                                s,
                                BoxBoundary,
                            )
                            G_fine = SolveRoutine(pointsBox)
                            # Compute std of qmc and expected value
                            QMC_std_next, QMC_Q =
                                ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                            if (SampleExponentCubature - 1 != 0)
                                if !cond(
                                    estimatedCubatureError,
                                    tolerance,
                                    isSearched,
                                    QMC_std_next,
                                )
                                    push!(estimatedCubatureErrorsInternals, QMC_std_next)
                                    push!(numberofSamples, numberOfPointsBox)
                                    push!(
                                        estimatedCubatureErrorsInternalsTimings,
                                        timingQMC +
                                        estimatedCubatureErrorsInternalsTimings[end],
                                    )
                                end


                            end
                        end


                    end




                    println(
                        "qmc error is ",
                        QMC_std,
                        " on box ",
                        -BoxBoundary,
                        " ",
                        BoxBoundary,
                        " exp val is ",
                        QMC_Q,
                        " exact sol is ",
                        analytical_sol(BoxBoundary, s, 2.6),
                        " samples is ",
                        numberOfPointsBox,
                    )
                    #SampleExponentCubature = SampleExponentCubature + 1
                    SampleExponentCubature = f(SampleExponentCubature, 1)

                end

                push!(SampleExponentCubatureArray, SampleExponentCubature - 1)

                println("------Finished Cubature error--------")

                push!(QMCResultsInternals, QMC_Q[1])

                # Routine to check evolution of truncation error
                if length(QMCResultsInternals) > 1
                    truncationError =
                        abs.(QMCResultsInternals[end-1] - QMCResultsInternals[end]) /
                        QMCResultsInternals[end]
                    println("truncation error is ", truncationError)

                    push!(estimatedTruncationErrorsInternals, truncationError)
                else
                    push!(estimatedTruncationErrorsInternals, -1)
                end

                if truncationError > tolerance / 2
                    estimatedCubatureErrorsInternals = [] # clear array when computing new truncation error
                    numberofSamples = []
                    lastTiming = estimatedCubatureErrorsInternalsTimings[end]
                    estimatedCubatureErrorsInternalsTimings = []
                    push!(estimatedCubatureErrorsInternalsTimings, lastTiming)

                end


                SampleExponentBox = SampleExponentBox + 1 # increase exponent for next box

            end
            SampleExponentBox = SampleExponentBox - 1   #if while loop of truncation error is satisfied reset the sample exponent back with one, avoid unecessary increasing the starting box for new tolerances 


            println("##############################################")


        end # end of @elapsed
        # SampleExponentBox = SampleExponentCubature
        #     end # end of if else


        #### computation of exact sol and updating arrays for printing

        DictOfEstimatedTruncationErrors[counter] = estimatedTruncationErrorsInternals
        DictOfEstimatedCubatureErrors[counter] = estimatedCubatureErrorsInternals
        DictOfSamples[counter] = numberofSamples




        exactSolOnBox = analytical_sol(BoxBoundary, s, 2.6)

        exactSol = analytical_sol(1000, s, 2.6)

        push!(
            exactCubatureErrors,
            abs(exactSolOnBox - QMCResultsInternals[end]) ./ QMCResultsInternals[end],
        )
        push!(exactTruncationErrors, abs(exactSol - exactSolOnBox) ./ exactSol)

        push!(estimatedCubatureErrors, estimatedCubatureErrorsInternals[end])
        push!(estimatedTruncationErrors, estimatedTruncationErrorsInternals[end])
        #push!(estimatedTime, t)

        if (length(estimatedTime)) > 0
            DictOfEstimatedCubatureErrorsTimings[counter] =
                estimatedCubatureErrorsInternalsTimings .+ estimatedTime[end]
            push!(estimatedTime, t + estimatedTime[end])
        else
            push!(estimatedTime, t)
            DictOfEstimatedCubatureErrorsTimings[counter] =
                estimatedCubatureErrorsInternalsTimings
        end

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
    loglog(estimatedTime, exactTruncationErrors, "*--k")
    loglog(estimatedTime, exactCubatureErrors, "*--r")

    loglog([3*10^-2, 3*10^-1], [10^-5, 10^-6], "--b")
    loglog([3*10^-2, 3*10^-1], [10^-5, 10^-7], "--m")
    loglog([3*10^-2, 3*10^-1], [10^-5, 10^-8], "--y")

    sz = 18

    grid(which = "both", ls = "-")
    legend(
        (
            "estimated truncation error",
            "estimated cubature error",
            "exact truncation error",
            "exact cubature error",
            "time^-1",
            "time^-2",
            "time^-3",
        ),
        fontsize = sz,
    )
    xlabel("time [sec]", fontsize = sz)
    ylabel("error", fontsize = sz)

    plt.xticks(fontsize = sz)
    plt.yticks(fontsize = sz)


    out=estimatedTime[1]
    DictOfEstimatedCubatureErrorsTimings[1][2]=out

    printy=Dict()
    printy=deepcopy(DictOfEstimatedCubatureErrors)
    printy[3][2]=printy[3][2]-0.00008



    for i = 1:length(tol)


        loglog(
            DictOfEstimatedCubatureErrorsTimings[i][2:end],
            DictOfEstimatedCubatureErrors[i][1:end],
            "-.or",
            alpha = 0.3,
            mec = "r",
        )

        for p = 1:1:length(DictOfSamples[i])
            num=log(DictOfSamples[i][p])/log(2)
            println(num)
            text(
                DictOfEstimatedCubatureErrorsTimings[i][2:end][p],
                printy[i][1:end][p],
                latexstring("\$2^{$num}\$"),
                fontsize = sz,rotation=0
            )
            println(latexstring("\2^{$num}"))
        end

    end

    println(DictOfEstimatedCubatureErrors[1][1:end])
    println(DictOfSamples[1])

    #    println(DictOfEstimatedCubatureErrorsTimings[1])
    #    println(DictOfEstimatedCubatureErrorsTimings[1][end])
    #    println(DictOfEstimatedCubatureErrorsTimings[2][end])
    #    println(estimatedTime)


    # Plot epoches i.e. the user requested tolerance
    for id = 1:length(estimatedTime)
        if (id == 1)
            loglog([0, estimatedTime[id]], [tol[id] / 2, tol[id] / 2], "--g")
            loglog([0, 0], [tol[id] * 10, tol[id] / 1000], "--g")
            loglog(
                [estimatedTime[id], estimatedTime[id]],
                [tol[id] * 10, tol[id] / 1000],
                "--g",
            )
        else
            loglog(
                [estimatedTime[id-1], estimatedTime[id]],
                [tol[id] / 2, tol[id] / 2],
                "--g",
            )
            loglog(
                [estimatedTime[id], estimatedTime[id]],
                [tol[id] * 10, tol[id] / 1000],
                "--g",
            )
        end

    end

    println(DictOfSamples)



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
    return QMC_std, QMC_Q[1]

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



function searchDirection(
    M::Int64,
    s::Int64,
    BoxBoundary::Float64,
    SampleExponentCubature::Int64,
    tolerance::Float64,
    QMCGenerator::Union{DigitalNet64,LatticeRule},
)
    """
    out = map(
        evalStd,
        [M, M, M],
        [s, s, s],
        [BoxBoundary, BoxBoundary, BoxBoundary],
        [SampleExponentCubature - 1, SampleExponentCubature, SampleExponentCubature + 1],
        [QMCGenerator, QMCGenerator, QMCGenerator],
    )
    """
    numberOfPointsBox = 2^(SampleExponentCubature+1)
    pointsBox = mapPoints(M, QMCGenerator, numberOfPointsBox, s, BoxBoundary)
    G_fine = SolveRoutine(pointsBox)

    out_1=ComputeQMCStdAndExp(G_fine[:,1:2^(SampleExponentCubature - 1),:],BoxBoundary,s,M)
    out0=ComputeQMCStdAndExp(G_fine[:,1:2^(SampleExponentCubature),:],BoxBoundary,s,M)
    out1=ComputeQMCStdAndExp(G_fine[:,1:2^(SampleExponentCubature+1),:],BoxBoundary,s,M)


    println(typeof(out_1))
    a=(out_1[1],out_1[2][1])
    b=(out0[1],out0[2][1])
    c=(out1[1],out1[2][1])

    #QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
    println(out_1[2][1])
    println(out0)
    println(out1)
    println("--- here ---")
    out=[a,b,c]
    println(out)

"""
    if out[2][1] < tolerance && out[1][1] > tolerance
        println("Choose ref value")
        return 1, out
    elseif out[2][1] < tolerance && out[1][1] < tolerance
        println("Backward search")
        return 2, out
    elseif out[2][1] > tolerance
        println("Forward search")
        return 3, out
    end
    """

    if out0[1] < tolerance && out_1[1] > tolerance
        println("Choose ref value")
        return 1, out
    elseif out_1[1] < tolerance && out0[1] < tolerance
        println("Backward search")
        return 2, out
    elseif out0[1] > tolerance
        println("Forward search")
        return 3, out
    end

end


function evalStd(
    M::Int64,
    s::Int64,
    BoxBoundary::Float64,
    SampleExponentCubature::Int64,
    QMCGenerator::Union{DigitalNet64,LatticeRule},
)


    numberOfPointsBox = 2^(SampleExponentCubature)
    pointsBox = mapPoints(M, QMCGenerator, numberOfPointsBox, s, BoxBoundary)


    G_fine = SolveRoutine(pointsBox)
    # Compute std of qmc and expected value
    QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)

    return QMC_std, QMC_Q[1]


end





main()
