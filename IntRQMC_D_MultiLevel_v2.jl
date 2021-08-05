using Random, Distributions
using DigitalNets
using LatticeRules
using PyPlot

#using MultilevelEstimators
#using LatticeRules: ShiftedLatticeRule

using SpecialFunctions: erf, erfinv

Φ⁻¹(x::T where {T<:Real}) = √2*erfinv(2*x-1)

Φ(x::T where {T<:Real}) = 1/2*(1+erf(x/√2))

function main()
    #GC.gc()

    #### Input parameters
    s = 2 # number of stochastic dimensions
    M = 4 # number of shifts
    N = 2  #2^N start number of samples
    b = -1:-0.5:-5.5
    RequestedTolerance_vec = 10 .^ b
    QMCGenerator = DigitalNet64(s)
    #QMCGenerator = LatticeRule(s)
    ####

    Data = RunSimulation(s, M, N, RequestedTolerance_vec, QMCGenerator)

    plotter(Data)

end # end of function main()


function RunSimulation(s::Int64, M::Int64, N::Int64, RequestedTolerance_vec::Vector, QMCGenerator::Random.AbstractRNG)

    if QMCGenerator isa AbstractDigitalNets
        fun_ShiftedQMCGenerator = (x) -> DigitalShiftedDigitalNets64(x)
        if Int64.(DigitalNets.reversebits((QMCGenerator.C)[1,1])) == 1
            println("Used point set is: Sobol_CS")
            QMCType = "Sobol"
        elseif Int64.(DigitalNets.reversebits((QMCGenerator.C)[1,1])) == 3
            println("Used point set is: Sobol_CS_2")
            QMCType = "Sobol3" ## FIXME?
        elseif Int64.(DigitalNets.reversebits((QMCGenerator.C)[1,1])) == 7
            println("Used point set is: Sobol_CS_3")
            QMCType = "Sobol3"
        end
    elseif QMCGenerator isa AbstractLatticeRule
        fun_ShiftedQMCGenerator = (x) -> ShiftedLatticeRule(x)
        println("Lattice")
        QMCType = "Lattice"
    end

    #Random.seed!(12345) # Setting the seed
    #funTransform = (x, y) -> transform(MultilevelEstimators.TruncatedNormal(0, 1, -y, y), x)

    f(x, gamma) = prod(1 .+ x .* gamma, dims=1)
    distrb = Distributions.Normal()

    RequestedTolerance_vec = vec(RequestedTolerance_vec)
    time_vec = zeros(length(RequestedTolerance_vec))
    err_vec = zeros(length(RequestedTolerance_vec))
    truncationError_vec = zeros(length(RequestedTolerance_vec))
    Real_truncationError_vec = zeros(length(RequestedTolerance_vec))
    Real_cubatureError_vec = zeros(length(RequestedTolerance_vec))

    cubatureError_vec = zeros(length(RequestedTolerance_vec))
    Nb_samples_vec = zeros(length(RequestedTolerance_vec))

    idx = 1
    while idx <= length(RequestedTolerance_vec)
        GC.gc()

        #Matrx_std=zeros(length(N),1)
        ct = 1 # FIXME: there is a second ct variable down there in the loop over Nshift

        level = 12
        sol = zeros(level+1)
        #sol_uncor=zeros(level+1)
        box_fine = 0
        box_coarse = 0
        Sample_Multiplier = 0
        l = 0
        truncationError = 10000
        cubatureError = 10000
        soll = 10000
        RequestedTolerance = RequestedTolerance_vec[idx]
        t = @elapsed begin
            percentagePtInSmallBox = []
            while l <= level
                p_f = 2^(N+l+Sample_Multiplier)
                p_c = 2^(N+l-1+Sample_Multiplier)

                p_f_FOR_BOX = 2^(N+l)
                p_c_FOR_BOX = 2^(N+l-1)

                println("---------")
                println("level ", l)
                println("Samples fine ", p_f)
                println("Samples coarse ", p_c)


                Matrix_f = zeros(s,p_f,M)
                Matrix_c = zeros(s,p_f,M)
                sizeShift = zeros(M)
                a = 1:1:s
                gamma = 1 ./ a .^2


                for Nshift = 1:M

                    ShiftedQMCGenerator = fun_ShiftedQMCGenerator(QMCGenerator)

                    ct = 0
                    for id = 1:length(ShiftedQMCGenerator[0:p_f_FOR_BOX]) - 1
                        box_fine = sqrt(2*2*log(p_f_FOR_BOX))
                        box_coarse = sqrt(2*2*log(p_c_FOR_BOX))

                        #Matrix_f[:,id,Nshift] = map.(funTransform, ShiftedQMCGenerator[id-1], box_fine)
                        Matrix_f[:,id,Nshift] = map.(x -> Φ⁻¹(Φ(-box_fine) + (Φ(box_fine) - Φ(-box_fine))*x), ShiftedQMCGenerator[id-1])

                        #cumsum = sum((Matrix_f[:,id,Nshift] .> box_coarse)  .| (Matrix_f[:,id,Nshift] .< -box_coarse))
                        #if cumsum > 0
                        if any((Matrix_f[:,id,Nshift] .> box_coarse)  .| (Matrix_f[:,id,Nshift] .< -box_coarse))
                            # point is not in the small box
                        else
                            # point is in the small box
                            ct = ct + 1
                            Matrix_c[:,ct,Nshift] = Matrix_f[:,id,Nshift]
                        end
                    end
                    sizeShift[Nshift] = ct
                end


                nshifts_effective = minimum(sizeShift)
                Matrix_c = Matrix_c[:,1:Int64(nshifts_effective),:] ## FIXME?
                pct = (1 - size(Matrix_c,2)/size(Matrix_f,2)) * 100
                println("Percentage of points in smaller box is ", pct, "%")
                println("box fine is ", -box_fine, " ", box_fine)
                println("box coarse is ", -box_coarse, " ", box_coarse)

                append!(percentagePtInSmallBox, pct)

                # Matrix_f is s-by-N-by-M --f--> 1-by-N-by-M
                G_fine = mean(f(Matrix_f, gamma), dims=2) # 1-by-1-by-M
                corrfactor_fine = (cdf(distrb, box_fine) - cdf(distrb, -box_fine))^s
                G_coarse = mean(f(Matrix_c, gamma), dims=2) # 1-by-1-by-M
                corrfactor_coarse = (cdf(distrb, box_coarse) - cdf(distrb, -box_coarse))^s
                diff = (G_fine .- G_coarse) * (corrfactor_fine - corrfactor_coarse)
                println("corrfactor big box ", corrfactor_fine)
                println("corrfactor small box ", corrfactor_coarse)
                println("corrfactor diff ", (corrfactor_fine - corrfactor_coarse))

                QMc_Q = mean(diff, dims=3) # estimate for truncation error over shifts; problem of number of shifts
                sol[l+1] = QMc_Q[1]

                #QMc_R = mean(G_fine * corrfactor_fine, dims=2) # = G_fine * corrfactor_fine
                QMc_R = G_fine * corrfactor_fine
                QMc_std = std(QMc_R) / sqrt(nshifts_effective) ## Not problem of number of shifts

                cubatureError = QMc_std
                truncationError = sol[l+1]

                println("Cubature error is ", cubatureError)
                println("Truncation error is ", sol[l+1])

                soll = mean(G_fine * corrfactor_fine, dims=3)[1]
                println("Sol is ", soll)
                println("abs error is ", abs(1-soll))

                if cubatureError > 0.5 * RequestedTolerance || truncationError == 0
                    Sample_Multiplier = Sample_Multiplier + 1
                else
                    if abs(truncationError) > 0.5 * RequestedTolerance && truncationError !=0
                        l = l + 1
                    else
                        l = level + 2
                        truncationError_vec[idx] = truncationError
                        cubatureError_vec[idx] = cubatureError
                        Real_truncationError_vec[idx] = 1 - corrfactor_fine
                        Real_cubatureError_vec[idx] = abs(corrfactor_fine - soll)
                    end
                end

                Nb_samples_vec[idx] = Nb_samples_vec[idx] + p_f

            end# end of for Nshift=1:M

        end # end of    while l<=level
        println("Runtime is ", t, " sec")
        println("******************************************************************************")

        time_vec[idx] = t
        err_vec[idx] = abs(1 - soll)
        idx = idx + 1
    end# end     while idx<=length(RequestedTolerance_vec)

    Data = Dict()
    Data[1] = RequestedTolerance_vec
    Data[2] = time_vec
    Data[3] = err_vec
    Data[4] = cubatureError_vec
    Data[5] = Real_cubatureError_vec
    Data[6] = truncationError_vec
    Data[7] = Real_truncationError_vec
    Data[8] = QMCType
    Data[9] = s
    Data[10] = M
    Data[11] = Nb_samples_vec

    return Data
end



function plotter(Data::Dict)

    RequestedTolerance_vec = Data[1]
    time_vec = Data[2]
    err_vec = Data[3]
    cubatureError_vec = Data[4]
    Real_cubatureError_vec = Data[5]
    truncationError_vec = Data[6]
    Real_truncationError_vec = Data[7]
    QMCType = Data[8]
    s = Data[9]
    M = Data[10]
    Nb_samples_vec = Data[11]

    figure()
    loglog(RequestedTolerance_vec, time_vec, "-*")
    str = string(QMCType, "\n", "s = ", s, "shift = ", M, " \n", "Runtime in function of requested tolerance")
    title(str)
    xlabel("RMSE")
    ylabel("run time sec")

    figure()
    loglog(RequestedTolerance_vec, err_vec, "-k*")
    println(time_vec)
    loglog(RequestedTolerance_vec, abs.(cubatureError_vec), "-r*")
    loglog(RequestedTolerance_vec, abs.(Real_cubatureError_vec), "--r*")
    loglog(RequestedTolerance_vec, abs.(truncationError_vec), "-g*")
    loglog(RequestedTolerance_vec, abs.(Real_truncationError_vec), "--g*")
    str = string(QMCType, "\n", "s = ", s, "shift = ", M, " \n", "Errors in function of requested RMSE")
    title(str)
    ylabel("Abs. error")
    xlabel("RMSE")
    legend(("abs. error on sol", "cubature error", "real cubature error", "abs. trunc. error", "abs. real trunc. error"))
    println(truncationError_vec)
    println(cubatureError_vec)

    figure()
    str = string(QMCType, "\n", "s = ", s, "shift = ", M, " \n", "Samples")
    title(str)
    loglog(RequestedTolerance_vec,Nb_samples_vec,"-*")
    ylabel("Nb. Samples")
    xlabel("Requested RMSE")
    println(Nb_samples_vec)
    println(RequestedTolerance_vec)

    #GC.gc()
end

main()
