using DigitalNets
using LatticeRules

using PyPlot

using Statistics: mean, std
using SpecialFunctions: erf, erfinv, gamma, gamma_inc
using StringLiterals

Φ⁻¹(x::T where {T<:Real}) = √2*erfinv(2*x-1)

Φ(x::T where {T<:Real}) = 1/2*(1+erf(x/√2))

#analytical_sol(a::Real,s::Int64) = ((gamma(4/5)/(2^(1/5))-gamma(4/5)*gamma_inc(4/5,a^2/2,0)[2]/(2^(1/5)))*2/sqrt(2*pi)+(Φ(a) - Φ(-a)))^s
analytical_sol(a::Real,s::Int64,params::Float64) =  params == 0.6 ? ((gamma(4/5)/(2^(1/5))-gamma(4/5)*gamma_inc(4/5,a^2/2,0)[2]/(2^(1/5)))*2/sqrt(2*pi)+(Φ(a) - Φ(-a)))^s : params == 2.6 ? ((gamma(9/5)*(2^(4/5))-gamma(9/5)*gamma_inc(9/5,a^2/2,0)[2]*(2^(4/5)))*2/sqrt(2*pi)+(Φ(a) - Φ(-a)))^s : 1



function main()

    #### Input parameters
    s = 3 # number of stochastic dimensions
    M = 16 # number of shifts
    alpha = 1
    N = 2 .^ collect(4:1:20)
    params = 0.6


    #generator = DigitalNet64InterlacedTwo(s)
    generator = LatticeRule(s)

    Data = RunSimulation(s, M, N, alpha, params, generator)

    #plotter(Data)

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
        RunSimulation(s, M, requestedTolerances, QMCGenerator)

        Run the algorithm for our test function in `s` dimensions with `M` shifts and a
            starting number of points `2^N0` for all the requested tolerances in
            `requestedTolerances` using the QMC point generator in `QMCGenerator`.
            """
            function RunSimulation(s::Int64, M::Int64, N::Vector, alpha::Int64, params::Float64, QMCGenerator::Union{DigitalNet64,LatticeRule})

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
                exactSol=0
                f(x) = prod(1 .+ abs.(x) .^ params, dims=1)


                # why do these two variables need to be declared here? can't we rewrite the logic?
                BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
                #truncationError = 10000
                #cubatureError = 10000
                soll = 10000
                idx=1
                for ell in N

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

                        # We first generate *all* the points for all the shifts...
                        BoxBoundary = sqrt(2*alpha*log(numberOfPointsBox))
                        boundsOfBoxes[idx] = BoxBoundary

                        for shiftId = 1:M
                            shiftedQMCGenerator = randomizedGenerator(QMCGenerator)

                            for id = 1:numberOfPointsBox

                                pointsBox[:,id,shiftId] = map.(x -> Φ⁻¹(Φ(-BoxBoundary) + (Φ(BoxBoundary) - Φ(-BoxBoundary))*x), shiftedQMCGenerator[id-1])

                            end
                        end

                        # pointsBox is s-by-N-by-M --f--> 1-by-N-by-M
                        G_fine = mean(f(pointsBox), dims=2) # 1-by-1-by-M
                        corrfactor_fine = (Φ(BoxBoundary) - Φ(-BoxBoundary))^s

                        QMC_R = abs.(G_fine) * (corrfactor_fine)


                        QMC_Q = mean(QMC_R, dims=3)
                        QMCResults[idx] = QMC_Q[1]


                        QMC_std = std(QMC_R) / sqrt(M)



                        cubature_error = QMC_std





                        exactSolOnBox = analytical_sol(BoxBoundary,s,params)
                        solutionsOnBox[idx] = exactSolOnBox
                        exactSol = analytical_sol(1000,s,params)
                        exactCubatureError = abs(exactSolOnBox - QMC_Q[1])./ QMC_Q[1]
                        exactTruncationError = abs(exactSol - exactSolOnBox)./exactSol
                        totError = abs(exactSol - QMC_Q[1])./exactSol


                        println("Levels needed ", ell)
                        println("samples needed on finest level ", ell)
                        println("box is ", BoxBoundary)
                        println("exact solution on given box is ", exactSolOnBox)
                        println("solution is ", QMC_Q[1])
                        println("exact solution  is ", exactSol)
                        println("Estimated rel Cubature error is ", cubature_error)
                        println("Exact rel Cubature error is ", exactCubatureError)
                        println("Estimated rel Truncation error is ", truncation_error)
                        println("Exact rel Truncation error is ", exactTruncationError)
                        println("total rel error is ", totError)




                    end # end of @elapsed

                    println("Runtime is ", t, " sec")
                    println("******************************************************************************")
                    estimatedTruncationErrors[idx] = truncation_error
                    estimatedCubatureErrors[idx] = cubature_error
                    timings[idx] = t
                    trueErrors[idx] = totError
                    trueTruncationErrors[idx] = exactTruncationError
                    trueCubatureErrors[idx] = exactCubatureError
                    maxLevelPerRequestedTol[idx] = ell
                    idx=idx+1

                end

                println(QMCType)
                println("alpha is " ,alpha)
                println("exact solution is ",exactSol)
                println("stochastic dimensions is equal to ",s)
                println("params is equal to ", params)
                println("number of shifts is " ,M)

                pr"""n   \%10s("a")   \%20s("Q")   \%25s("tot rel err")   \%15s("Iab") \%20s("trunc rel err") \%20s("box cub rel err") \n"""
                idx=1
                for i in N
                    pr"\%-3d(N[idx]) \%16.8f(boundsOfBoxes[idx])       \%16.16f(QMCResults[idx])       \%.5e(trueErrors[idx])        \%16.16f(solutionsOnBox[idx])       \%.5e(trueTruncationErrors[idx])       \%.5e(trueCubatureErrors[idx])\n"

                    idx=idx+1


                end


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
                loglog(requestedTolerances, requestedTolerances.^2, "--k+")
                loglog(requestedTolerances, requestedTolerances, "--r+")


                str = string(QMCType, "\n", "s = ", s, " shift = ", M, " \n", "Errors in function of requested RMSE")
                title(str)
                ylabel("Abs. error")
                xlabel("RMSE")
                legend(("true error", "est cubature error", "true cubature error", "est trunc error", "true trunc error"))
                println(estimatedTruncationErrors)
                println(estimatedCubatureErrors)
                println(trueTruncationErrors)
                println(trueCubatureErrors)

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
