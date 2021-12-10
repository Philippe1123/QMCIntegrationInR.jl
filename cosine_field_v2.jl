
using PyPlot
using FiniteElementDiffusion
using DelimitedFiles
using SpecialFunctions: erf, erfinv, gamma, gamma_inc


Φ⁻¹(x::T where {T<:Real}) = √2 * erfinv(2 * x - 1)
Φ⁻¹(x::T where {T<:Real},σ::T where {T<:Real}) = √2 * erfinv(2 * x - 1) * σ

Φ(x::T where {T<:Real}) = 1 / 2 * (1 + erf(x / √2))
Φ(x::T where {T<:Real},σ::T where {T<:Real}) = 1 / 2 * (1 + erf(x /(σ * √2)))
cdfNorm(x::T where {T<:Real}) = 1 / (sqrt(2 * pi)) .* exp.(-(x .^ 2) ./ 2)



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
    
    NumberOfElementsPlus1=1000
    pts = collect(0:1/(NumberOfElementsPlus1):1)



    power = 2


  #  figure("a")
  #  figure("b")
   # figure("field")
    a=-5
    b=5
        utot=0
        maxpt=2^13
        s=1

        u_vec=zeros(maxpt,1)
        sample_vec=zeros(maxpt,s)
        Field_vec=zeros(maxpt,1)
        figure("field")
        for op=1:maxpt
            samplesPoints = ((b .- a).*rand(s,1) .+ a)
  #          samplesPoints = randn(s,1) 

            n = 1
            Field=zeros(length(pts),1)
            for l in samplesPoints
                Field=Field .+ l./n.^power .* cos.(2*pi*n*pts)
                n=n+1
            end
            plot(Field)
            sleep(2)
            Field =  (Field)
            Field_vec[op]=Field[5]

 



            sample_vec[op,1]=samplesPoints[1]
            sample_vec[op,2]=samplesPoints[2]

      




        end


       figure("projectie 1")
       plot(sample_vec[:,1],vec(Field_vec),"*")
        figure("projectie 2")
        plot(sample_vec[:,2],vec(Field_vec),"*")
        figure("field 3d")
        surf(sample_vec[:,1],sample_vec[:,2],vec(Field_vec))





       

end




main()