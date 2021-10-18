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

Φ⁻¹(x::T where {T<:Real}) = √2 * erfinv(2 * x - 1)

Φ(x::T where {T<:Real}) = 1 / 2 * (1 + erf(x / √2))





analytical_sol(a::Real, s::Int64, sigma::Real) =
    (
        (
            gamma((1 + sigma) / 2) -
            gamma((1 + sigma) / 2) * gamma_inc((1 + sigma) / 2, a^2 / 2, 0)[2]
        ) * 2^(sigma / 2) / sqrt(pi) + erf(a / sqrt(2))
    )^s

function main()

    #### Input parameters
    M = 64 # number of shifts
    alpha = 3
    N_lattice = 2 .^ collect(4:1:20)
    N_net = 2 .^ collect(4:1:20)


    #  generator = DigitalNet64InterlacedTwo(s)
    #generator = DigitalNet64InterlacedThree(s)

    #generator = LatticeRule(s)

    #generator = DigitalNet64(s)

    
    # dim = 1
    s = 1

    
    Data = RunSimulation(
        s,
        M,
        N_lattice,
        1,
        0.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "1")
#############################################################################    
    
    Data = RunSimulation(s, M, N_net, 1, 0.6, DigitalNet64(s))
    writeOut(Data, "2")

#############################################################################    


    Data = RunSimulation(
        s,
        M,
        N_lattice,
        2,
        1.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "3")
#############################################################################    


    Data = RunSimulation(s, M, N_net, 2, 1.6, DigitalNet64InterlacedTwo(s))
    writeOut(Data, "4")
#############################################################################    

    Data = RunSimulation(
        s,
        M,
        N_lattice,
        3,
        2.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "5")
#############################################################################    

    Data = RunSimulation(s, M, N_net, 3, 2.6, DigitalNet64(readdlm("nx_s5_alpha3_m32.col",UInt64),s,2^30))
    writeOut(Data, "6")   

#############################################################################    
    
    # dim = 2
    s = 2
    Data = RunSimulation(
        s,
        M,
        N_lattice,
        1,
        0.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "7")
#############################################################################    

    Data = RunSimulation(s, M, N_net, 1, 0.6, DigitalNet64(s))
    writeOut(Data, "8")
#############################################################################    

    Data = RunSimulation(
        s,
        M,
        N_lattice,
        2,
        1.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "9")
#############################################################################    

    Data = RunSimulation(s, M, N_net, 2, 1.6, DigitalNet64InterlacedTwo(s))
    writeOut(Data, "10")
#############################################################################    

    Data = RunSimulation(
        s,
        M,
        N_lattice,
        3,
        2.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "11")
#############################################################################    

    Data = RunSimulation(s, M, N_net, 3, 2.6, DigitalNet64InterlacedThree(s))
    writeOut(Data, "12")
#############################################################################    
      

    # dim = 3
    s = 3
    Data = RunSimulation(
        s,
        M,
        N_lattice,
        1,
        0.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "13")
#############################################################################    

    Data = RunSimulation(s, M, N_net, 1, 0.6, DigitalNet64(s))
    writeOut(Data, "14")

#############################################################################    

    Data = RunSimulation(
        s,
        M,
        N_lattice,
        2,
        1.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "15")
#############################################################################    

    Data = RunSimulation(s, M, N_net, 2, 1.6, DigitalNet64InterlacedTwo(s))
    writeOut(Data, "16")
#############################################################################    

    Data = RunSimulation(
        s,
        M,
        N_lattice,
        3,
        2.6,
        LatticeRule(vec(UInt32.(readdlm("exew_base2_m20_a3_HKKN.txt"))), s),
    )
    writeOut(Data, "17")
#############################################################################    

    Data = RunSimulation(s, M, N_net, 3, 2.6, DigitalNet64InterlacedThree(s))
    writeOut(Data, "18")
#############################################################################    


end # end of function main()


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
    shiftAverages = zeros(M, length(N))



    #f(x, params) = prod(
    #    (1 .+ abs.(x) .^ params) .* 1 / (sqrt(2 * pi)) .* exp.(-(x .^ 2) ./ 2),
    #    dims = 1,
    #)

    f(x, params) = prod(
        (1 .+ abs.(x) .^ params) .* 1 / (sqrt(2 * pi)) .* exp.(-(x .^ 2) ./ 2),
        dims = 1,
    )


    # why do these two variables need to be declared here? can't we rewrite the logic?
    BoxBoundary = 0 # our large box is [-largeBox, largeBox]^d
    #truncationError = 10000
    #cubatureError = 10000
    soll = 10000
    idx = 1
    for ell in N
        println("Currently running sample number ", ell, " exponent ", log.(ell) ./ log(2))
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





           




            # pointsBox is s-by-N-by-M --f--> 1-by-N-by-M
            G_fine = mean(f(pointsBox, params), dims = 2) # 1-by-1-by-M
            pointsBox = 0
            GC.gc()




            QMC_R = abs.(G_fine) * (BoxBoundary * 2)^s
            G_fine = 0
            GC.gc()
            shiftAverages[:, idx] = vec(reshape(QMC_R, M, 1, 1))

       



            QMC_Q = mean(QMC_R, dims = 3)
   

            QMCResults[idx] = QMC_Q[1]


            QMC_std = std(QMC_R) / sqrt(M)
            QMC_R = 0
            GC.gc()


            cubature_error = QMC_std





            exactSolOnBox = analytical_sol(BoxBoundary, s, params)
            solutionsOnBox[idx] = exactSolOnBox
            exactSol = analytical_sol(1000, s, params)
            exactCubatureError = abs(exactSolOnBox - QMC_Q[1]) ./ QMC_Q[1]
            exactTruncationError = abs(exactSol - exactSolOnBox) ./ exactSol
            totError = abs(exactSol - QMC_Q[1]) ./ exactSol


            # println("Levels needed ", ell)
            # println("samples needed on finest level ", ell)
            # println("box is ", BoxBoundary)
            # println("exact solution on given box is ", exactSolOnBox)
            # println("solution is ", QMC_Q[1])
            # println("exact solution  is ", exactSol)
            # println("Estimated rel Cubature error is ", cubature_error)
            # println("Exact rel Cubature error is ", exactCubatureError)
            # println("Estimated rel Truncation error is ", truncation_error)
            # println("Exact rel Truncation error is ", exactTruncationError)
            # println("total rel error is ", totError)




        end # end of @elapsed

        #println("Runtime is ", t, " sec")
        #println(
        #    "******************************************************************************",
        #)
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
        ", alpha = ",
        Data[14],
    )
    title(str)
    ylabel("Rel. error")
    xlabel("N")
    legend((
        "true error",
        "true cubature error",
        "true truncation error",
        "N^-1",
        "N^-2",
        "N^-3",
    ))

    savefig(string(str, ".png"))
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")
    println("    ")

    figure()
    plt.boxplot(
        Data[15],
        positions = Int64.(log.(Data[11]) ./ log(2)),# Each column/cell is one box
        notch = true, # Notched center
        whis = 0.75, # Whisker length as a percent of inner quartile range
        widths = 0.25, # Width of boxes
        vert = true, # Horizontal boxes
        sym = "rs",
    ) # Symbol color and shape (rs = red square)
    #  plt.yscale("log")
    #  plt.xscale("log")

    title(str)
    ylabel("expected value")
    xlabel("log 2 of number of samples")
    savefig(string(str, "boxplot.png"))

    figure()
    plt.boxplot(
        Data[15][:, (end-4):end],
        positions = Int64.(log.(Data[11][(end-4):end]) ./ log(2)),# Each column/cell is one box
        notch = true, # Notched center
        whis = 0.75, # Whisker length as a percent of inner quartile range
        widths = 0.25, # Width of boxes
        vert = true, # Horizontal boxes
        sym = "rs",
    ) # Symbol color and shape (rs = red square)
    #  plt.yscale("log")
    #  plt.xscale("log")

    title(str)
    ylabel("expected value")
    xlabel("log 2 of number of samples")
    savefig(string(str, "boxplot_zoom.png"))


end

function writeOut(Data::Dict, str::String)

    currentFolder = @__DIR__
    Path = string(currentFolder, "/", str, ".jld2")
    println(Path)
    save(Path, "data", Data)


end

main()
