using DigitalNets
using LatticeRules

#using PyPlot

using Statistics
using SpecialFunctions: erf, erfinv, gamma, gamma_inc
using StringLiterals
using PrettyTables
using Random
using DelimitedFiles
using JLD2
using FileIO
using PyPlot
using FiniteElementDiffusion
using GaussianRandomFields



NormalPdf(x) = 1 / sqrt(2 * pi) * exp(-1 / 2 * x^2)



function main()


    s = 3
    M = 8


    tol = 10.0 .^ (-1:-1:-8)

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







    BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
    SampleExponentBox = 2
    counter = 1

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

        if length(estimatedTruncationErrors) > 0 && estimatedTruncationErrors[end] < tolerance / 2 &&estimatedCubatureErrors[end] < tolerance / 2
           

                estimatedTruncationErrorsInternals = DictOfEstimatedTruncationErrors[counter-1] 
                estimatedCubatureErrorsInternals = DictOfEstimatedCubatureErrors[counter-1] 
                estimatedCubatureErrorsInternalsTimings = DictOfEstimatedCubatureErrorsTimings[counter-1]
                
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

                        # compute all the points
                        pointsBox = mapPoints(
                            M,
                            QMCGenerator,
                            numberOfPointsBox,
                            s,
                            BoxBoundary,
                        )

                        # Solve the problem
                        G_fine = SolveRoutine(pointsBox, s)
                        QMC_std, QMC_Q = ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
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
                    println(
                        "qmc error is ",
                        QMC_std,
                        " on box ",
                        -BoxBoundary,
                        " ",
                        BoxBoundary,
                        " exp val is ",
                        QMC_Q,
                    )




                    ####################### Adjust Cubature error
                    estimatedCubatureError = QMC_std
                    SampleExponentCubature = SampleExponentBox    # start with exponent used for box
                    while estimatedCubatureError > tolerance / 2
                        timingQMC = @elapsed begin

                            numberOfPointsBox = 2^(SampleExponentCubature)
                            pointsBox = mapPoints(
                                M,
                                QMCGenerator,
                                numberOfPointsBox,
                                s,
                                BoxBoundary,
                            )
                            G_fine = SolveRoutine(pointsBox, s)
                            # Compute std of qmc and expected value
                            QMC_std, QMC_Q =
                                ComputeQMCStdAndExp(G_fine, BoxBoundary, s, M)
                            estimatedCubatureError = QMC_std
                            push!(estimatedCubatureErrorsInternals, QMC_std)

                        end
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
                            QMC_Q,
                        )
                        SampleExponentCubature = SampleExponentCubature + 1

                    end
                    println("------Finished Cubature error--------")

                    push!(QMCResultsInternals, QMC_Q[1])

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

    grid(which = "both", ls = "-")
    legend((
        "estimatedTruncationErrors",
        "estimatedCubatureErrors",
        "time^-1",
        "time^-2",
        "time^-3",
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


function SolveRoutine(pointsBox::Array, s::Int64)



    MaterialParam = Dict()
    QuadPoints = 3

    #Order 1
    Elements=Int64.(readdlm(joinpath(locationOfMesh,"2D/Structured/Quad/Elements_1_16.txt")))
    Elements = Elements[:, 5:end]
    Nodes=readdlm(joinpath(locationOfMesh,"2D/Structured/Quad/Nodes_1_16.txt"))
    Nodes1=Nodes[:,2:3]#only retain xy component

    Center=compute_centers(Nodes1,Elements)
    matField = GaussianRandomFields.Matern(0.3,2.0,σ=1.0,p=2)
    cov = CovarianceFunction(2,matField)
    grf =  GaussianRandomField(cov,KarhunenLoeve(s),Center,Elements[:,1:3],quad=GaussLegendre())



    ElemType="TwoD_Quad_Order1"
    NumberOfElements = size(Elements, 1)

    G_fine = zeros(1, size(pointsBox, 2), size(pointsBox, 3))


    # Define random field Gaussian random field
   
    for j = 1:size(pointsBox, 2) #loop over samples
        for k = 1:size(pointsBox, 3) #loop over shifts
            samplesPoints = pointsBox[:, j, k]

            Field = GaussianRandomFields.sample(grf, xi = samplesPoints)

            Field = 0.1 .+ exp.(Field)
            # fem routine
            for id = 1:NumberOfElements
                MaterialParam[id] = Field[id]
            end
            solverparam = (
                elemtype = ElemType,
                Qpt = QuadPoints,
                Nelem = NumberOfElements,
                Order = parse(Int, ElemType[end]),
            )
            u1 = solver2D.main(Nodes1, Elements, MaterialParam, solverparam)
            # fem routine end
            #select mid point

            u1 = u1[Int64(ceil(length(u1) / 2))] * prod(NormalPdf.((samplesPoints)))
            G_fine[1, j, k] = u1
        end
    end
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

function compute_centers(p, t)
    d = size(p, 2)
    vec_t = vec(t)
    size_t = size(t)

    pts = Array{Float64}(undef, size(t, 1), d)
    @inbounds for i = 1:d
        x = reshape(p[vec_t, i], size_t)
        mean!(view(pts, :, i), x)
    end
    pts
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
