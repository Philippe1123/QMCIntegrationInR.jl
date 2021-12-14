
using PyPlot
using FiniteElementDiffusion
using DelimitedFiles
using SpecialFunctions: erf, erfinv, gamma, gamma_inc


Φ⁻¹(x::T where {T<:Real}) = √2 * erfinv(2 * x - 1)
Φ⁻¹(x::T where {T<:Real},σ::T where {T<:Real}) = √2 * erfinv(2 * x - 1) * σ

Φ(x::T where {T<:Real}) = 1 / 2 * (1 + erf(x / √2))
Φ(x::T where {T<:Real},σ::T where {T<:Real}) = 1 / 2 * (1 + erf(x /(σ * √2)))
NormalPdf(x) = 1/sqrt(2*pi) *exp(-1/2*x^2)

logNormalPdf(x) = x <= 0 ? 0 : (1 / (x*sqrt(2 * pi)) .* exp.(-(log(x) ).^ 2 ./ 2))


function main()
    MaterialParam = Dict()
    QuadPoints=6

    #Order 1
    Elements=Int64.(readdlm(joinpath(locationOfMesh,"1D/Elements_1_5.txt")))
    Elements=Elements[:,5:end]
    Nodes=readdlm(joinpath(locationOfMesh,"1D/Nodes_1_5.txt"))
    Nodes1=Nodes[:,2]#only retain x component
    ElemType="OneD_Order1"
    NumberOfElements=size(Elements,1)
    
    NumberOfElementsPlus1=10
    pts = collect(0:1/(NumberOfElementsPlus1):1)



    power = 2


  #  figure("a")
  #  figure("b")
   # figure("field")
    a=-6
    b=6
        utot=0
        maxpt=2^20
        s=1
     #   figure("u1")
      #  figure("u1_mult")
        u_vec=zeros(maxpt,1)
        sample_vec=zeros(maxpt,s)
        Field_vec=zeros(maxpt,1)
        for op=1:maxpt
            samplesPoints = (((b .- a).*rand(s,1) .+ a))
            n = 1
            Field=zeros(length(pts),1)
            for l in samplesPoints
                Field=Field .+ l./n.^power .* sin.(pi*pts)
                n=n+1
            end
            FieldIn =  0.1.+exp.(Field[2:2:end])
            Field_vec[op]=Field[5]

            for id=1:NumberOfElements MaterialParam[id]=FieldIn[id] end
            solverparam=(elemtype =ElemType, Qpt=QuadPoints, Nelem=NumberOfElements, Order=parse(Int,ElemType[end]))
            u1_full=solver1D.main(Nodes1,Elements,MaterialParam,solverparam)
            """
            println()
            println(u1_full[Int64(floor(length(u1_full)/2))] )
            println(logNormalPdf(Field[5]))
            println(logNormalPdf(samplesPoints[1]))


            println(u1_full[Int64(floor(length(u1_full)/2))] * logNormalPdf(samplesPoints[1]))
            println("---")
            sleep(2.)
            """
            
            u1=u1_full[Int64(floor(length(u1_full)/2))] * prod(NormalPdf.(samplesPoints))
            """
            figure("u1")
            plot(u1_full,"-*")
            figure("u1_mult")
            plot(u1_full .*  logNormalPdf(samplesPoints[1]) + u1_full .* logNormalPdf(samplesPoints[2]),"-*")

            sleep(2)
            """

            u_vec[op]=u1
            sample_vec[op,1]=(samplesPoints[1])
            #sample_vec[op,2]=(samplesPoints[2])

            utot=u1+utot
            
         #  sleep(4)
         #  figure("Field").clear()
         #  figure("dist").clear()
         #  figure("u1_mod").clear()
         #  figure("u1").clear()



        end
        
        figure("projectie 1")
        plot(sample_vec[:,1],vec(u_vec),"*")
       #surf(sample_vec[:,1],sample_vec[:,2],vec(u_vec))
        #figure("projectie 2")
        #plot(sample_vec[:,2],vec(u_vec),"*")
       # figure("3d")
       # plot(sample_vec[:,1],vec(u_vec))
       # figure("field 3d")
      #  plot(sample_vec[:,1],vec(Field_vec))s
        
        println(utot/maxpt*(b-a)^s)
       # println(utot/maxpt)




        u_vec=zeros(maxpt,1)
        sample_vec=zeros(maxpt,1)
        BoxBoundary=100
        utot=0
        out = -1
        for op=1:maxpt
            samplesPoints = rand(s,1)

            out=map.(
                x -> Φ⁻¹(
                    Φ(-BoxBoundary,1) + (Φ(BoxBoundary,1) - Φ(-BoxBoundary,1)) * x
                ),
                samplesPoints,
            )
            n = 1
            Field=zeros(length(pts),1)
            for l in out
                Field=Field .+ l./n.^power .* sin.(pi*pts)
                n=n+1
            end
            Field =  0.1.+exp.(Field)
            FieldIn =  Field[2:2:end]

            for id=1:NumberOfElements MaterialParam[id]=FieldIn[id] end
            solverparam=(elemtype =ElemType, Qpt=QuadPoints, Nelem=NumberOfElements, Order=parse(Int,ElemType[end]))
            u1=solver1D.main(Nodes1,Elements,MaterialParam,solverparam)
          #  figure("u1")
          #  plot(pts[2:2:end-2],u1,"-*")
            u1=u1[2] 
           # figure("Field")
          #  plot(pts,Field,"-*")

            utot=u1+utot
            u_vec[op]=u1
            sample_vec[op]=samplesPoints[1]
          #  sleep(1)
          

        end
        println(utot/maxpt)
        figure("projection with mapping")
        plot(vec(sample_vec),vec(u_vec),"*")



end




main()