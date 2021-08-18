using DigitalNets
using LatticeRules

using PyPlot

using Statistics: mean, std
using SpecialFunctions: erf, erfinv

Φ⁻¹(x::T where {T<:Real}) = √2*erfinv(2*x-1)

Φ(x::T where {T<:Real}) = 1/2*(1+erf(x/√2))

function main()

    #### Input parameters
    s = 2 # number of stochastic dimensions
    M = 16 # number of shifts
    N0 = 2  #2^N0 start number of samples
    b = -1:-0.5:-5
    requestedTolerances = 10 .^ b

#    generator = DigitalNet64InterlacedTwo(s)
    generator = LatticeRule(s)

    Data = RunSimulation(s, M, N0, requestedTolerances, generator)

    plotter(Data)

end # end of function main()


## Quick and dirty solution for the "isa" problem (needs to be fixed in a more decent matter)
"""
Create a randomized generator with a random shift or random digital shift for the passed in QMC generator.
    """
    randomizedGenerator(latticeGenerator::LatticeRule) = ShiftedLatticeRule(latticeGenerator)
    randomizedGenerator(digitalnetGenerator::DigitalNet64) = DigitalShiftedDigitalNets64(digitalnetGenerator)

    ## Quick and dirty solution for labeling the point set (needs to be fixed in a more decent matter)
    """
    String label for the current generator, either Lattice or Sobol1, Sobol2, or Sobol3.
        """
        labelThisGenerator(latticeGenerator::LatticeRule) = "Lattice"
        labelThisGenerator(digitalnetGenerator::DigitalNet64) = "Sobol$(Int32(round(sqrt(reversebits(digitalnetGenerator.C[1])))))"


        """
        RunSimulation(s, M, N0, requestedTolerances, QMCGenerator)

        Run the algorithm for our test function in `s` dimensions with `M` shifts and a
            starting number of points `2^N0` for all the requested tolerances in
            `requestedTolerances` using the QMC point generator in `QMCGenerator`.
            """
            function RunSimulation(s::Int64, M::Int64, N0::Int64, requestedTolerances::Vector, QMCGenerator::Union{DigitalNet64,LatticeRule})

                QMCType = labelThisGenerator(QMCGenerator)

                # This is our current integrand function, it should be an argument to this function... to be fixed
                a = 1:1:s
                gamma = 1 ./ a .^2
                f(x) = prod(1 .+ x .* gamma, dims=1)

                timings = zeros(length(requestedTolerances))
                trueErrors = zeros(length(requestedTolerances))
                estimatedTruncationErrors = zeros(length(requestedTolerances))
                estimatedCubatureErrors = zeros(length(requestedTolerances))
                trueTruncationErrors = zeros(length(requestedTolerances))
                trueCubatureErrors = zeros(length(requestedTolerances))
                nbOfSamples = zeros(length(requestedTolerances))

                maxLevel = 20

                for (idx, requestedTolerance) in enumerate(requestedTolerances)

                    sol = zeros(maxLevel+1)
                    # why do these two variables need to be declared here? can't we rewrite the logic?
                    largeBoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
                    sampleMultiplierLog2 = 0
                    ell = 0
                    #truncationError = 10000
                    #cubatureError = 10000
                    soll = 10000

                    t = @elapsed begin

                        solutions_per_level = zeros(1,maxLevel)
                        cubature_error = 1000
                        truncation_error = 1000
                        QMc_Q = 0


                        while ell <= maxLevel
                            cubature_error = 1000
                            truncation_error = 1000


                            pLargeBox = 2^(ell)
                            N = 0;
                            QMc_Q = 0


                            while cubature_error > 0.5 * requestedTolerance

                                numberOfPointsLargeBox = 2^(N+ell+sampleMultiplierLog2)
                                pointsLargeBox = zeros(s, numberOfPointsLargeBox, M)

                                println("---------")
                                println("level ", ell)
                                println("p for large box ", pLargeBox)
                                println("Samples  ", N)

                                # We first generate *all* the points for all the shifts...
                                # This does not seem like a very good idea.
                                for shiftId = 1:M
                                    shiftedQMCGenerator = randomizedGenerator(QMCGenerator)
                                    largeBoxBoundary = sqrt(2*2*log(pLargeBox))

                                    for id = 1:numberOfPointsLargeBox # Note: this said pLargeBox instead of numberOfPointsLargeBox

                                        pointsLargeBox[:,id,shiftId] = map.(x -> Φ⁻¹(Φ(-largeBoxBoundary) + (Φ(largeBoxBoundary) - Φ(-largeBoxBoundary))*x), shiftedQMCGenerator[id-1])

                                    end
                                end

                                # pointsLargeBox is s-by-N-by-M --f--> 1-by-N-by-M
                                G_fine = mean(f(pointsLargeBox), dims=2) # 1-by-1-by-M
                                corrfactor_fine = (Φ(largeBoxBoundary) - Φ(-largeBoxBoundary))^s

                                vector_of_solutions = abs.(G_fine) * (corrfactor_fine)
                                println("corrfactor big box ", corrfactor_fine)


                                QMc_Q = mean(vector_of_solutions, dims=3)

                                #QMc_R = mean(G_fine * corrfactor_fine, dims=2) # = G_fine * corrfactor_fine
                                QMc_R = G_fine * corrfactor_fine
                                QMc_std = std(QMc_R) / sqrt(M)

                                cubature_error = QMc_std
                                N=N+1

                            end
                            println("Cubature error is ", cubature_error)

                            if ell == 0
                                solutions_per_level[ell+1] = QMc_Q[1]
                                ell = ell + 1
                            else
                                solutions_per_level[ell+1] = QMc_Q[1]
                                truncation_error = abs(solutions_per_level[ell+1] - solutions_per_level[ell])
                                println("Truncation error is ", truncation_error)

                                if truncation_error > 0.5 * requestedTolerance
                                    ell = ell + 1
                                else

                                    ell = maxLevel + 1

                                end

                            end

                        end # end of while l <= level

                    end # end of @elapsed

                    println("Runtime is ", t, " sec")
                    println("******************************************************************************")
                    estimatedTruncationErrors[idx] = truncation_error
                    estimatedCubatureErrors[idx] = cubature_error
                    timings[idx] = t
                    trueErrors[idx] = abs(1 - QMc_Q[1])
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
                str = string(QMCType, "\n", "s = ", s, " shift = ", M, " \n", "Runtime in function of requested tolerance")
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
                str = string(QMCType, "\n", "s = ", s, " shift = ", M, " \n", "Errors in function of requested RMSE")
                title(str)
                ylabel("Abs. error")
                xlabel("RMSE")
                legend(("true error", "est cubature error", "true cubature error", "est trunc error", "true trunc error"))
                println(estimatedTruncationErrors)
                println(estimatedCubatureErrors)

                figure()
                str = string(QMCType, "\n", "s = ", s, " shift = ", M, " \n", "Samples")
                title(str)
                loglog(requestedTolerances,nbOfSamples,"-*")
                ylabel("Nb. Samples")
                xlabel("Requested RMSE")
                println(nbOfSamples)
                println(requestedTolerances)

            end

            main()
