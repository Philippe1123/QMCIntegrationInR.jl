using Random, Distributions
using DigitalNets
using LatticeRules
using MultilevelEstimators
using PyPlot









function main()
    #    Random.seed!(12345) # Setting the seed
    #d = MultilevelEstimators.Normal()
    #d = MultilevelEstimators.TruncatedNormal(0.0,1.0,-2.0,2.0)


    DictOfMat = Dict()
    DictOfQMC_pt = Dict()
    funTransform = (x, y) -> transform(MultilevelEstimators.TruncatedNormal(0, 1, -y, y), x)

    id = 2

    N_dim = 2 .^ id
    N_small = 2 .^ (id - 1)


    BigBoxInterval = sqrt(2 * 2 * log(N_dim))
    SmallBoxInterval = sqrt(2 * 2 * log(N_small))

    N = N_dim * 2^2

    println("Sample points big box is ", N)
    println("Sample points small box is ", N_small)

    println("Big box is ", BigBoxInterval)
    println("Small box is ", SmallBoxInterval)

    Matrx_std = zeros(length(N), 1)
    ct = -2000
    s = 2
    Num_Shift = 1

    for p in N
        Mat = zeros(s, p, Num_Shift)
        Mat_full = zeros(s, p, Num_Shift)
        QMC_pt = zeros(s, p, Num_Shift)
        QMC_pt_intern = zeros(s, p, Num_Shift)

        Lattice = DigitalNet64_2(s)
        for Nshift = 1:Num_Shift
            shiftLat = DigitalShiftedDigitalNets64(Lattice)
            ct = 1
            for mp = 1:p
                #    Mat[:,ct,1]=map.(funTransform,Lattice[mp-1],BigBoxInterval)
                Mat_full[:, mp, Nshift] = map.(funTransform, shiftLat[mp-1], BigBoxInterval)
                QMC_pt[:, mp, Nshift] = shiftLat[mp-1]
                cumsum = sum(
                    (Mat_full[:, mp, Nshift] .> SmallBoxInterval) .|
                    (Mat_full[:, mp, Nshift] .< -SmallBoxInterval),
                )

                #println((Mat_full[:,mp,Nshift]))
                #println((Mat_full[:,mp,Nshift] .> SmallBoxInterval))
                #println((Mat_full[:,mp,Nshift] .< -SmallBoxInterval))
                #println((Mat_full[:,mp,Nshift] .> SmallBoxInterval)  .| (Mat_full[:,mp,Nshift] .< -SmallBoxInterval))

                if (cumsum > 0)
                else
                    Mat[:, ct, Nshift] = Mat_full[:, mp, Nshift]
                    QMC_pt_intern[:, ct, Nshift] = QMC_pt[:, mp, Nshift]
                    #cat(2,Mat,Mat_full[:,mp,1])
                    ct = ct + 1
                end


            end
            DictOfMat[Nshift] = Mat[:, 1:ct-1, Nshift]
            DictOfQMC_pt[Nshift] = QMC_pt_intern[:, 1:ct-1, Nshift]
        end
        #Mat=Mat[:,1:ct-1,:]
        #QMC_pt_intern=QMC_pt_intern[:,1:ct-1,:]
        #Mat=Mat[:,1:ct-1,:]
        #QMC_pt_intern=QMC_pt_intern[:,1:ct-1,:]
        println(ct)
        println(p)
        println("percentage lost ", (1 - (ct - 1) / p) * 100, "%")
        println("-------------")
        println("-------------")
        println("-------------")

        #println(size(Mat,2))
        #println(size(Mat_full,2))
        #println("Percentage of points in bigger box is ", (1-size(Mat,2)/size(Mat_full,2))*100)

        figure()
        plot(Mat_full[1, :, 1], Mat_full[2, :, 1], "r*")
        plot(DictOfMat[1][1, :], DictOfMat[1][2, :], "g*")
        plot(
            [SmallBoxInterval, SmallBoxInterval],
            [-SmallBoxInterval, SmallBoxInterval],
            "k",
        )

        plot(
            [-SmallBoxInterval, SmallBoxInterval],
            [SmallBoxInterval, SmallBoxInterval],
            "k",
        )
        plot(
            [-SmallBoxInterval, SmallBoxInterval],
            [-SmallBoxInterval, -SmallBoxInterval],
            "k",
        )

        plot(
            [-SmallBoxInterval, -SmallBoxInterval],
            [-SmallBoxInterval, SmallBoxInterval],
            "k",
        )
        plot(
            [SmallBoxInterval, -SmallBoxInterval],
            [SmallBoxInterval, SmallBoxInterval],
            "k",
        )

        strtitle = string(
            "Dimensioning number of samples is ",
            N_dim,
            "\n",
            " Actual number of samples  is ",
            N,
            "\n",
            "big box is [- ",
            round(BigBoxInterval, digits = 2),
            " ",
            round(BigBoxInterval, digits = 2),
            "] ",
            "small box is [- ",
            round(SmallBoxInterval, digits = 2),
            " ",
            round(SmallBoxInterval, digits = 2),
            "]",
        )
        title(strtitle)
        #legend(("Samples included in difference of big and small box","Samples only included in small box"))

        figure()
        plot(QMC_pt[1, :, 1], QMC_pt[2, :, 1], "r*")
        plot(DictOfQMC_pt[1][1, :], DictOfQMC_pt[1][2, :], "g*")
        title(strtitle)
        #legend(("Samples included in difference of big and small box","Samples only included in small box"))

    end





end


main()
