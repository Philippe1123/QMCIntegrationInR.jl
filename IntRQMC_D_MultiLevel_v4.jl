using DigitalNets
using LatticeRules

using PyPlot

using Statistics: mean, std
using SpecialFunctions: erf, erfinv, gamma, gamma_inc
using StringLiterals
using PrettyTables

Φ⁻¹(x::T where {T<:Real}) = √2 * erfinv(2 * x - 1)

Φ(x::T where {T<:Real}) = 1 / 2 * (1 + erf(x / √2))



#analytical_sol(a::Real,s::Int64) = ((gamma(4/5)/(2^(1/5))-gamma(4/5)*gamma_inc(4/5,a^2/2,0)[2]/(2^(1/5)))*2/sqrt(2*pi)+(Φ(a) - Φ(-a)))^s
analytical_sol(a::Real, s::Int64, params::Float64) =
    params == 0.6 ?
    (
        (
            gamma(4 / 5) / (2^(1 / 5)) -
            gamma(4 / 5) * gamma_inc(4 / 5, a^2 / 2, 0)[2] / (2^(1 / 5))
        ) * 2 / sqrt(2 * pi) + (Φ(a) - Φ(-a))
    )^s :
    params == 2.6 ?
    (
        (
            gamma(9 / 5) * (2^(4 / 5)) -
            gamma(9 / 5) * gamma_inc(9 / 5, a^2 / 2, 0)[2] * (2^(4 / 5))
        ) * 2 / sqrt(2 * pi) + (Φ(a) - Φ(-a))
    )^s : 1



function main()

    #### Input parameters
    s = 1 # number of stochastic dimensions
    M = 16 # number of shifts
    alpha = 2
    N = 2 .^ collect(4:1:17)
    params = 0.6
    correctionFactor = true


    generator = DigitalNet64InterlacedTwo(s)
    #generator = LatticeRule(s)

    Data = RunSimulation(s, M, N, alpha, params, generator, correctionFactor)

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
    f(x) = prod(1 .+ abs.(x) .^ params, dims = 1)


    # why do these two variables need to be declared here? can't we rewrite the logic?
    BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
    #truncationError = 10000
    #cubatureError = 10000
    soll = 10000
    idx = 1
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
            BoxBoundary = sqrt(2 * alpha * log(numberOfPointsBox))
            boundsOfBoxes[idx] = BoxBoundary

            for shiftId = 1:M
                shiftedQMCGenerator = randomizedGenerator(QMCGenerator)

                for id = 1:numberOfPointsBox
                    pointsBox[:, id, shiftId] =
                        map.(
                            x -> Φ⁻¹(
                                Φ(-BoxBoundary) + (Φ(BoxBoundary) - Φ(-BoxBoundary)) * x,
                            ),
                            shiftedQMCGenerator[id-1],
                        )


                end
            end

            # pointsBox is s-by-N-by-M --f--> 1-by-N-by-M
            G_fine = mean(f(pointsBox), dims = 2) # 1-by-1-by-M

            corrfactor_fine =
                correctionFactor == true ? (Φ(BoxBoundary) - Φ(-BoxBoundary))^s : 1

            QMC_R = abs.(G_fine) * (corrfactor_fine)


            QMC_Q = mean(QMC_R, dims = 3)
            QMCResults[idx] = QMC_Q[1]


            QMC_std = std(QMC_R) / sqrt(M)



            cubature_error = QMC_std





            exactSolOnBox = analytical_sol(BoxBoundary, s, params)
            solutionsOnBox[idx] = exactSolOnBox
            exactSol = analytical_sol(1000, s, params)
            exactCubatureError = abs(exactSolOnBox - QMC_Q[1]) ./ QMC_Q[1]
            exactTruncationError = abs(exactSol - exactSolOnBox) ./ exactSol
            totError = abs(exactSol - QMC_Q[1]) ./ exactSol


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
        println(
            "******************************************************************************",
        )
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
        "n",
        "a",
        "Q",
        "tot rel error",
        "Iab",
        "trunc rel error (exact)",
        "box cub rel error (exact)",
        "corr factor ",
    ])

    formatters = ft_printf(
        ["%-3d", "%16.8f", "%16.16f", "%.5e", "%16.16f", "%.5e", "%.5e", "%16.8f"],
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
    Data[12] = correctionFactor
    Data[13] = params
    Data[14] = alpha

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
        ", corr factor ? = ",
        Data[12],
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

main()
