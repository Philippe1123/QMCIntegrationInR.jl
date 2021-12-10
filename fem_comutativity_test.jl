
using PyPlot
using FiniteElementDiffusion
using DelimitedFiles
using SpecialFunctions: erf, erfinv, gamma, gamma_inc


Φ⁻¹(x::T where {T<:Real}) = √2 * erfinv(2 * x - 1)
Φ⁻¹(x::T where {T<:Real},σ::T where {T<:Real}) = √2 * erfinv(2 * x - 1) * σ

Φ(x::T where {T<:Real}) = 1 / 2 * (1 + erf(x / √2))
Φ(x::T where {T<:Real},σ::T where {T<:Real}) = 1 / 2 * (1 + erf(x /(σ * √2)))
cdfNorm(x::T where {T<:Real}) = 1 / (sqrt(2 * pi)) .* exp.(-(x .^ 2) ./ 2)

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
    a=-2
    b=2
        utot=0
        maxpt=1
        s=2
        figure("u1")
        figure("u1_mult")
        u_vec=zeros(maxpt,1)
        sample_vec=zeros(maxpt,s)
        Field_vec=zeros(maxpt,1)
        for op=1:maxpt
            samplesPoints = ((b .- a).*rand(s,1) .+ a)
            n = 1
            Field=zeros(length(pts),1)
            Field_part=Dict()
            for l in samplesPoints
                Field=Field .+ l./n.^power .* sin.(pi*pts)
                Field_part[n] = l./n.^power .* sin.(pi*pts)
                n=n+1
            end

            println(Field)
            println(Field_part[1])
            println(Field_part[2])
            println(Field_part[1]+Field_part[2])

            #full field
            FieldIn =  exp.(Field[2:2:end])
            Field_vec[op]=Field[5]
            for id=1:NumberOfElements MaterialParam[id]=FieldIn[id] end
            solverparam=(elemtype =ElemType, Qpt=QuadPoints, Nelem=NumberOfElements, Order=parse(Int,ElemType[end]))
            u1_full=solver1D.main(Nodes1,Elements,MaterialParam,solverparam)

            #part field 1
            Field = Field_part[1]
            FieldIn =  exp.(Field[2:2:end])
            Field_vec[op]=Field[5]
            for id=1:NumberOfElements MaterialParam[id]=FieldIn[id] end
            solverparam=(elemtype =ElemType, Qpt=QuadPoints, Nelem=NumberOfElements, Order=parse(Int,ElemType[end]))
            u2_full=solver1D.main(Nodes1,Elements,MaterialParam,solverparam)



            #part field 2

            Field = Field_part[2]
            FieldIn =  exp.(Field[2:2:end])
            Field_vec[op]=Field[5]
            for id=1:NumberOfElements MaterialParam[id]=FieldIn[id] end
            solverparam=(elemtype =ElemType, Qpt=QuadPoints, Nelem=NumberOfElements, Order=parse(Int,ElemType[end]))
            u3_full=solver1D.main(Nodes1,Elements,MaterialParam,solverparam)



            println(u1_full)
            println(u2_full .+ u3_full)

            figure()
            plot(u1_full,"-r")
            plot(u2_full .+ u3_full,"-b")
            plot(u2_full,"--b")
            plot(u3_full,"-.b")



        end
        
        



end




main()