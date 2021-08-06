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
    b = -1:-0.5:-5.5
    requestedTolerances = 10 .^ b

    #generator = DigitalNet64(s)
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

    maxLevel = 12

    for (idx, requestedTolerance) in enumerate(requestedTolerances)

        sol = zeros(maxLevel+1)
        # why do these two variables need to be declared here? can't we rewrite the logic?
        largeBoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
        smallBoxBoundary = 0 # our small box is [-smallBox, smallBox]^d
        sampleMultiplierLog2 = 0
        ell = 0
        #truncationError = 10000
        #cubatureError = 10000
        soll = 10000

        t = @elapsed begin

            percentagePtInSmallBox = []

            while ell <= maxLevel

                pLargeBox = 2^(N0+ell)
                pSmallBox = 2^(N0+ell-1)

                numberOfPointsLargeBox = 2^(N0+ell+sampleMultiplierLog2)
                #numberOfPointsSmallBox = 2^(N0+ell-1+sampleMultiplierLog2)

                println("---------")
                println("level ", ell)
                println("p for large box ", pLargeBox)
                println("p for small box ", pSmallBox)
                println("Samples for large box ", numberOfPointsLargeBox)
                #println("Samples coarse ", numberOfPointsSmallBox)

                pointsLargeBox = zeros(s, numberOfPointsLargeBox, M)
                pointsSmallBox = zeros(s, numberOfPointsLargeBox, M)
                nbOfPointsInSmallBoxPerShift = zeros(M)

                # We first generate *all* the points for all the shifts...
                # This does not seem like a very good idea.
                for shiftId = 1:M

                    shiftedQMCGenerator = randomizedGenerator(QMCGenerator)

                    largeBoxBoundary = sqrt(2*2*log(pLargeBox))
                    smallBoxBoundary = sqrt(2*2*log(pSmallBox))

                    ct = 0 # count the points which fall into the small box
                    for id = 1:numberOfPointsLargeBox # Note: this said pLargeBox instead of numberOfPointsLargeBox

                        pointsLargeBox[:,id,shiftId] = map.(x -> Φ⁻¹(Φ(-largeBoxBoundary) + (Φ(largeBoxBoundary) - Φ(-largeBoxBoundary))*x), shiftedQMCGenerator[id-1])

                        if any((pointsLargeBox[:,id,shiftId] .> smallBoxBoundary)  .| (pointsLargeBox[:,id,shiftId] .< -smallBoxBoundary))
                            # point is not in the small box
                        else
                            # point is in the small box
                            ct = ct + 1
                            pointsSmallBox[:,ct,shiftId] = pointsLargeBox[:,id,shiftId]
                        end
                    end

                    nbOfPointsInSmallBoxPerShift[shiftId] = ct

                end

                nshiftsEffective = minimum(nbOfPointsInSmallBoxPerShift)
                pointsSmallBox = pointsSmallBox[:,1:Int64(nshiftsEffective),:] ## FIXME?
                pct = (1 - size(pointsSmallBox,2)/size(pointsLargeBox,2)) * 100
                println("Percentage of points in smal box is ", pct, "%")
                println("Large box is ", -largeBoxBoundary, " ", largeBoxBoundary)
                println("Small box is ", -smallBoxBoundary, " ", smallBoxBoundary)

                append!(percentagePtInSmallBox, pct)

                # pointsLargeBox is s-by-N-by-M --f--> 1-by-N-by-M
                G_fine = mean(f(pointsLargeBox), dims=2) # 1-by-1-by-M
                corrfactor_fine = (Φ(largeBoxBoundary) - Φ(-largeBoxBoundary))^s
                G_coarse = mean(f(pointsSmallBox), dims=2) # 1-by-1-by-M
                corrfactor_coarse = (Φ(smallBoxBoundary) - Φ(-smallBoxBoundary))^s
                diff = abs.(G_fine .- G_coarse) * (corrfactor_fine - corrfactor_coarse)
                println("corrfactor big box ", corrfactor_fine)
                println("corrfactor small box ", corrfactor_coarse)
                println("corrfactor diff ", (corrfactor_fine - corrfactor_coarse))

                QMc_Q = mean(diff, dims=3) # estimate for truncation error over shifts; problem of number of shifts
                sol[ell+1] = QMc_Q[1]

                #QMc_R = mean(G_fine * corrfactor_fine, dims=2) # = G_fine * corrfactor_fine
                QMc_R = G_fine * corrfactor_fine
                QMc_std = std(QMc_R) / sqrt(nshiftsEffective) ## Not problem of number of shifts

                cubatureError = QMc_std
                truncationError = sol[ell+1]

                println("Cubature error is ", cubatureError)
                println("Truncation error is ", sol[ell+1])

                soll = mean(G_fine * corrfactor_fine, dims=3)[1]
                println("Sol is ", soll)
                println("abs error is ", abs(1-soll))

                if cubatureError > 0.5 * requestedTolerance || truncationError == 0
                    sampleMultiplierLog2 = sampleMultiplierLog2 + 1
                else
                    if abs(truncationError) > 0.5 * requestedTolerance && truncationError != 0
                        ell = ell + 1
                    else
                        ell = maxLevel + 2 # such that we stop the while loop
                        estimatedTruncationErrors[idx] = truncationError
                        estimatedCubatureErrors[idx] = cubatureError
                        trueTruncationErrors[idx] = 1 - corrfactor_fine
                        trueCubatureErrors[idx] = abs(corrfactor_fine - soll)
                    end
                end

                nbOfSamples[idx] = nbOfSamples[idx] + numberOfPointsLargeBox

            end # end of while l <= level

        end # end of @elapsed

        println("sampleMultiplierLog2 is ", sampleMultiplierLog2)
        println("Runtime is ", t, " sec")
        println("******************************************************************************")

        timings[idx] = t
        trueErrors[idx] = abs(1 - soll)
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
