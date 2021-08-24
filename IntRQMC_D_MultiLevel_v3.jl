using DigitalNets
using LatticeRules

using PyPlot

using Statistics: mean, std
using SpecialFunctions: erf, erfinv, gamma, gamma_inc

Φ⁻¹(x::T where {T<:Real}) = √2 * erfinv(2 * x - 1)

Φ(x::T where {T<:Real}) = 1 / 2 * (1 + erf(x / √2))

#analytical_sol(a::Real,s::Int64) = ((gamma(4/5)/(2^(1/5))-gamma(4/5)*gamma_inc(4/5,a^2/2,0)[2]/(2^(1/5)))*2/sqrt(2*pi)+(Φ(a) - Φ(-a)))^s
analytical_sol(a::Real, s::Int64) =
    (
        (
            gamma(9 / 5) * (2^(4 / 5)) -
            gamma(9 / 5) * gamma_inc(9 / 5, a^2 / 2, 0)[2] * (2^(4 / 5))
        ) * 2 / sqrt(2 * pi) + (Φ(a) - Φ(-a))
    )^s

function main()

    #### Input parameters
    s = 1 # number of stochastic dimensions
    M = 16 # number of shifts
    N0 = 2  #2^N0 start number of samples
    b = -1:-0.5:-5
    requestedTolerances = 10 .^ b

    generator = DigitalNet64InterlacedTwo(s)
    #    generator = LatticeRule(s)

    Data = RunSimulation(s, M, N0, requestedTolerances, generator)

    plotter(Data)

end # end of function main()


## Quick and dirty solution for the "isa" problem (needs to be fixed in a more decent matter)
"""
Create a randomized generator with a random shift or random digital shift for the passed in QMC generator.
    """
randomizedGenerator(latticeGenerator::LatticeRule) = ShiftedLatticeRule(latticeGenerator)
randomizedGenerator(digitalnetGenerator::DigitalNet64) =
    DigitalShiftedDigitalNets64(digitalnetGenerator)

## Quick and dirty solution for labeling the point set (needs to be fixed in a more decent matter)
"""
String label for the current generator, either Lattice or Sobol1, Sobol2, or Sobol3.
    """
labelThisGenerator(latticeGenerator::LatticeRule) = "Lattice"
labelThisGenerator(digitalnetGenerator::DigitalNet64) =
    "Sobol$(Int32(round(sqrt(reversebits(digitalnetGenerator.C[1])))))"


"""
RunSimulation(s, M, N0, requestedTolerances, QMCGenerator)

Run the algorithm for our test function in `s` dimensions with `M` shifts and a
    starting number of points `2^N0` for all the requested tolerances in
    `requestedTolerances` using the QMC point generator in `QMCGenerator`.
    """
function RunSimulation(
    s::Int64,
    M::Int64,
    N0::Int64,
    requestedTolerances::Vector,
    QMCGenerator::Union{DigitalNet64,LatticeRule},
)

    QMCType = labelThisGenerator(QMCGenerator)

    # This is our current integrand function, it should be an argument to this function... to be fixed
    a = 1:1:s
    gamma = 1 ./ a .^ 2
    #                f(x) = prod(1 .+ x .* gamma, dims=1)
    #        f(x) = prod(1 .+ abs.(x) .^ 0.6, dims=1)
    f(x) = prod(1 .+ abs.(x) .^ 2.6, dims = 1)

    timings = zeros(length(requestedTolerances))
    trueErrors = zeros(length(requestedTolerances))
    estimatedTruncationErrors = zeros(length(requestedTolerances))
    estimatedCubatureErrors = zeros(length(requestedTolerances))
    trueTruncationErrors = zeros(length(requestedTolerances))
    trueCubatureErrors = zeros(length(requestedTolerances))
    nbOfSamples = zeros(length(requestedTolerances))
    maxLevelPerRequestedTol = Int64.(zeros(length(requestedTolerances)))

    maxLevel = 30

    for (idx, requestedTolerance) in enumerate(requestedTolerances)
        sol = zeros(maxLevel + 1)
        # why do these two variables need to be declared here? can't we rewrite the logic?
        BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
        ell = 0
        #truncationError = 10000
        #cubatureError = 10000
        soll = 10000

        t = @elapsed begin

            solutions_per_level = zeros(1, maxLevel + 1)
            cubature_error = 1000 # initalisation
            truncation_error = 1000 # initalisation
            QMC_Q = 0
            corrfactor_fine = 0
            BoxBoundary = 0
            exactTruncationError = 1000 # initalisation
            exactCubatureError = 1000 # initalisation
            totError = 1000

            if idx == 1
                ell = 0
            else

                ell = maxLevelPerRequestedTol[idx-1]
            end

            while ell <= maxLevel
                cubature_error = 1000  # reset
                truncation_error = 1000 # reset
                exactTruncationError = 1000 # reset
                exactCubatureError = 1000 # reset
                totError = 1000


                pBox = 2^(ell)
                N = 0
                QMC_Q = 0


                while cubature_error > 0.5 * requestedTolerance

                    numberOfPointsBox = 2^(N)
                    pointsBox = zeros(s, numberOfPointsBox, M)

                    # We first generate *all* the points for all the shifts...
                    # This does not seem like a very good idea.
                    for shiftId = 1:M
                        shiftedQMCGenerator = randomizedGenerator(QMCGenerator)
                        BoxBoundary = sqrt(2 * 2 * log(pBox))

                        for id = 1:numberOfPointsBox
                            pointsBox[:, id, shiftId] =
                                map.(
                                    x -> Φ⁻¹(
                                        Φ(-BoxBoundary) +
                                        (Φ(BoxBoundary) - Φ(-BoxBoundary)) * x,
                                    ),
                                    shiftedQMCGenerator[id-1],
                                )

                        end
                    end

                    # pointsBox is s-by-N-by-M --f--> 1-by-N-by-M
                    G_fine = mean(f(pointsBox), dims = 2) # 1-by-1-by-M
                    corrfactor_fine = (Φ(BoxBoundary) - Φ(-BoxBoundary))^s

                    QMC_R = abs.(G_fine) * (corrfactor_fine)


                    QMC_Q = mean(QMC_R, dims = 3)

                    #QMC_R = mean(G_fine * corrfactor_fine, dims=2) # = G_fine * corrfactor_fine
                    #QMC_R = G_fine * corrfactor_fine
                    QMC_std = std(QMC_R) / sqrt(M)



                    cubature_error = QMC_std
                    N = N + 1

                end

                if ell == 0
                    solutions_per_level[ell+1] = QMC_Q[1]
                    ell = ell + 1
                else
                    solutions_per_level[ell+1] = QMC_Q[1]
                    truncation_error =
                        abs(solutions_per_level[ell+1] - solutions_per_level[ell])


                    if truncation_error > 0.5 * requestedTolerance
                        ell = ell + 1
                    else

                        exactSolOnBox = analytical_sol(BoxBoundary, s)
                        exactSol = analytical_sol(1000, s)
                        exactCubatureError = abs(exactSolOnBox - QMC_Q[1])
                        exactTruncationError = abs(exactSol - exactSolOnBox)
                        totError = abs(exactSol - QMC_Q[1])


                        println("For tolerance ", requestedTolerance)
                        println("Levels needed ", ell)
                        println("samples needed on finest level ", 2^N)
                        println("box is ", BoxBoundary)
                        println("exact solution on given box is ", exactSolOnBox)
                        println("solution is ", QMC_Q[1])
                        println("exact solution  is ", exactSol)
                        println("Estimated Cubature error is ", cubature_error)
                        println("Exact Cubature error is ", exactCubatureError)
                        println("Estimated Truncation error is ", truncation_error)
                        println("Exact Truncation error is ", exactTruncationError)
                        println("total  error is ", totError)

                        break


                    end

                end

            end # end of while l <= level

        end # end of @elapsed

        println("Runtime is ", t, " sec")
        println(
            "******************************************************************************",
        )
        estimatedTruncationErrors[idx] = truncation_error
        estimatedCubatureErrors[idx] = cubature_error
        timings[idx] = t
        trueErrors[idx] = totError
        trueTruncationErrors[idx] = exactTruncationError
        trueCubatureErrors[idx] = exactCubatureError
        maxLevelPerRequestedTol[idx] = ell
    end # loop over requestedTolerances



    Data = Dict()
    Data[1] = requestedTolerances
    Data[2] = timings
    Data[3] = trueErrors
    Data[4] = estimatedCubatureErrors
    Data[5] = trueCubatureErrors
    Data[6] = estimatedTruncationErrors
    Data[7] = trueTruncationErrors
    Data[8] = QMCType
    Data[9] = s
    Data[10] = M
    Data[11] = nbOfSamples

    return Data
end



function plotter(Data::Dict)

    requestedTolerances = Data[1]
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
    loglog(requestedTolerances, timings, "-*")
    str = string(
        QMCType,
        "\n",
        "s = ",
        s,
        " shift = ",
        M,
        " \n",
        "Runtime in function of requested tolerance",
    )
    title(str)
    xlabel("RMSE")
    ylabel("run time sec")

    figure()
    loglog(requestedTolerances, trueErrors, "-k*")
    println(timings)
    loglog(requestedTolerances, abs.(estimatedCubatureErrors), "-r*")
    loglog(requestedTolerances, abs.(trueCubatureErrors), "--r*")
    loglog(requestedTolerances, abs.(estimatedTruncationErrors), "-g*")
    loglog(requestedTolerances, abs.(trueTruncationErrors), "--g*")
    loglog(requestedTolerances, requestedTolerances .^ 2, "--k+")
    loglog(requestedTolerances, requestedTolerances, "--r+")


    str = string(
        QMCType,
        "\n",
        "s = ",
        s,
        " shift = ",
        M,
        " \n",
        "Errors in function of requested RMSE",
    )
    title(str)
    ylabel("Abs. error")
    xlabel("RMSE")
    legend((
        "true error",
        "est cubature error",
        "true cubature error",
        "est trunc error",
        "true trunc error",
    ))
    println(estimatedTruncationErrors)
    println(estimatedCubatureErrors)
    println(trueTruncationErrors)
    println(trueCubatureErrors)

    figure()
    str = string(QMCType, "\n", "s = ", s, " shift = ", M, " \n", "Samples")
    title(str)
    loglog(requestedTolerances, nbOfSamples, "-*")
    ylabel("Nb. Samples")
    xlabel("Requested RMSE")
    println(nbOfSamples)
    println(requestedTolerances)

end

main()
