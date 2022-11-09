using DigitalNets
using LatticeRules

using Random

using PyPlot
using SpecialFunctions
using DelimitedFiles

Φ⁻¹(x::T where {T<:Real}) = √2 * erfinv(2 * x - 1)
Φ⁻¹(x::T where {T<:Real},σ::T where {T<:Real}) = √2 * erfinv(2 * x - 1) * σ


function main()



   lat =  LatticeRule(2)
   pt=lat[0:2^10-1]
    figure()
    ptmat=Array{String}(undef, 2^10,3)


   for id=1:size(pt)[1]
    plot(pt[id][1],pt[id][2],"k*")
    ptmat[id,1]=string(pt[id][1])
    ptmat[id,2]=string(pt[id][2])
    ptmat[id,3]="\\\\"

   end

   stp = open(string(@__DIR__,"/pt.txt"), "w")
   writedlm(stp, ptmat)
   close(stp)

   println(ptmat)



   uniform_pt= mapPointsUniform(1,lat,2^10-1,2,1000.)


   println(size(uniform_pt))
   println(uniform_pt[1,:,1])
   println(uniform_pt[2,:,1])
   figure()
   plot(uniform_pt[1,:,1],uniform_pt[2,:,1],"k*")

   ptmatU=Array{String}(undef, 2^10-1,3)
   for id=1:size(uniform_pt,2)
    ptmatU[id,1]=string(uniform_pt[1,id,1])
    ptmatU[id,2]=string(uniform_pt[2,id,1])
    ptmatU[id,3]="\\\\"
   end
   println(ptmatU)


   stp = open(string(@__DIR__,"/ptmatU.txt"), "w")
   writedlm(stp, ptmatU)
   close(stp)

   gpt = mapPointsGaussian(1,lat,2^20-1,2,1)
   println(size(gpt))
    figure()
   plot(gpt[1,:,1],gpt[2,:,1],"k*")

   ptmatG1=Array{String}(undef, 2^10,3)
   for id=1:size(gpt,2)
    ptmatG1[id,1]=string(gpt[1,id,1])
    ptmatG1[id,2]=string(gpt[2,id,1])
    ptmatG1[id,3]="\\\\"
   end
   println(ptmatG1)

   stp = open(string(@__DIR__,"/ptmatG1.txt"), "w")
   writedlm(stp, ptmatG1)
   close(stp)



   gpt = mapPointsGaussian(1,lat,2^20-1,2,2)
   println(size(gpt))
    figure()
   plot(gpt[1,:,1],gpt[2,:,1],"k*")


   gpt = mapPointsGaussian(1,lat,2^10-1,2,3)
   println(size(gpt))
    figure()
   plot(gpt[1,:,1],gpt[2,:,1],"k*")

end
randomizedGenerator(latticeGenerator::LatticeRule) = ShiftedLatticeRule(latticeGenerator)
function randomizedGenerator(digitalnetGenerator::DigitalNet64)
    DigitalShiftedDigitalNets64(digitalnetGenerator)
end


function mapPointsUniform(
    M::Int64,
    QMCGenerator::Union{DigitalNet64,LatticeRule},
    numberOfPointsBox::Int64,
    s::Int64,
    BoxBoundary::Float64,
)

    Random.seed!(1234)
    pointsBox = zeros(s, numberOfPointsBox, M)
    for shiftId = 1:M
#        shiftedQMCGenerator = randomizedGenerator(QMCGenerator)
        shiftedQMCGenerator = (QMCGenerator)

        for id = 1:numberOfPointsBox
            pointsBox[:, id, shiftId] =
                map.(
                    x -> -BoxBoundary + (BoxBoundary - (-BoxBoundary)) * x,
                    shiftedQMCGenerator[id-1],
                )
        end
    end
    return pointsBox
end

function mapPointsGaussian(
    M::Int64,
    QMCGenerator::Union{DigitalNet64,LatticeRule},
    numberOfPointsBox::Int64,
    s::Int64,
    a::Int64
)
pointsBox = zeros(s, numberOfPointsBox, M)


Random.seed!(1234)
            for shiftId = 1:M
                shiftedQMCGenerator = (QMCGenerator)

                for id = 1:numberOfPointsBox
                    pointsBox[:, id, shiftId] =
                        map.(
                            x -> Φ⁻¹( x,1.
                            )*a,
                            shiftedQMCGenerator[id-1],
                        )

                end
            end
            return pointsBox

        end


main()
