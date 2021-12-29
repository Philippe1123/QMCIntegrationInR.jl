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
    M = 8


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







    BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
    SampleExponentBox = 2


    for tolerance in tol
        #SampleExponentBox = 2
        println("##############################################")
        println("Currently running tolerance number ", tolerance)
        println("")
        truncationError = 10 * tolerance # initalisation
        QMCResultsInternals = [] # reintilized for each tolerance, only used internally
        estimatedCubatureErrorsInternals = [] # reintilized for each tolerance, only used internally
        estimatedTruncationErrorsInternals = [] # reintilized for each tolerance, only used internally
        samplesInternals = [] # reintilized for each tolerance, only used internally
        boundsOfBoxesInternals = [] # reintilized for each tolerance, only used internally
        t = @elapsed begin
            ########################################## Adjust truncation error error
            while truncationError > tolerance / 2
                BoxBoundary = 0
                numberOfPointsBox = 2^(SampleExponentBox)
                push!(samplesInternals, numberOfPointsBox)

                BoxBoundary = sqrt(2 * alpha * log(numberOfPointsBox))
                push!(boundsOfBoxesInternals, BoxBoundary)


                pointsBox = mapPoints(M,QMCGenerator,numberOfPointsBox,s,BoxBoundary)

                # Solve the problem



                G_fine = SolveRoutine(pointsBox)
                QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                push!(estimatedCubatureErrorsInternals, QMC_std)



                ####################### Adjust Cubature error
                estimatedCubatureError = QMC_std
                SampleExponentCubature = SampleExponentBox
                while estimatedCubatureError > tolerance / 2
                    numberOfPointsBox = 2^(SampleExponentCubature)
                    pointsBox = mapPoints(M,QMCGenerator,numberOfPointsBox,s,BoxBoundary)
                    G_fine = SolveRoutine(pointsBox)
                    # Compute std of qmc and expected value
                    QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                    estimatedCubatureError = QMC_std
                    push!(estimatedCubatureErrorsInternals, QMC_std)
                    println("qmc error is ",QMC_std, " on box ", -BoxBoundary," ",BoxBoundary)
                    SampleExponentCubature = SampleExponentCubature + 1

                end
                println("------Finished Cubature error--------")

                push!(QMCResultsInternals, QMC_Q[1])

                # Routine to check evolution of truncation error
                if length(QMCResultsInternals) > 1
                    truncationError =
                        abs.(QMCResultsInternals[end-1] - QMCResultsInternals[end]) /
                        QMCResultsInternals[end]
                    println("truncation error is ",truncationError)

                    push!(estimatedTruncationErrorsInternals, truncationError)
                else
                    push!(estimatedTruncationErrorsInternals, -1)
                end

                SampleExponentBox = SampleExponentBox + 1

            end
                SampleExponentBox = SampleExponentBox - 1


            println("##############################################")


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


function mapPoints(M::Int64,QMCGenerator::Union{DigitalNet64,LatticeRule},numberOfPointsBox::Int64,s::Int64,BoxBoundary::Float64)

    Random.seed!(1234)
    pointsBox = zeros(s, numberOfPointsBox, M)
    for shiftId = 1:M
        shiftedQMCGenerator = randomizedGenerator(QMCGenerator)
        for id = 1:numberOfPointsBox
            pointsBox[:, id, shiftId] =
                map.(
                    x ->
                        -BoxBoundary + (BoxBoundary - (-BoxBoundary)) * x,
                    shiftedQMCGenerator[id-1],
                )
        end
    end
    return pointsBox
end






main()
