using DigitalNets
using GaussianRandomFields 
using PyPlot
using MultilevelEstimators
using LatticeRules
using FiniteElementDiffusion
using DelimitedFiles
using Statistics
using SpecialFunctions


function main()






MaterialParam=Dict()
DiffusionCoefficient=1.0
QuadPoints=6

#Order 1
Elements=Int64.(readdlm(joinpath(locationOfMesh,"1D/Elements_1_5.txt")))
Elements=Elements[:,5:end]
Nodes=readdlm(joinpath(locationOfMesh,"1D/Nodes_1_5.txt"))
Nodes1=Nodes[:,2]#only retain x component
ElemType="OneD_Order1"
NumberOfElements=size(Elements,1)


########## in Rd
eval_points=0+0.0001:0.0001:10
output=zeros(length(eval_points),1)
ct=1
for j in eval_points
for id=1:NumberOfElements MaterialParam[id]=j end
solverparam=(elemtype =ElemType, Qpt=QuadPoints, Nelem=NumberOfElements, Order=parse(Int,ElemType[end]))
u1=solver1D.main(Nodes1,Elements,MaterialParam,solverparam)
u1=u1[Int64(floor(length(u1)/2))] * 1 / (j*sqrt(2 * pi)) .* exp.(-(log(j) .^ 2) ./ 2)
output[ct]=u1*(10-0.0001)
ct=ct+1
end
println(mean(output))
figure()
plot(eval_points,output,"*")


########## in [0:1]
eval_points=0+0.00001:0.00001:1-0.00001
output=zeros(length(eval_points),1)
ct=1
a=1
for j in eval_points
    y = a*exp(erfinv(2*j-1)*sqrt(2))
    dy = a*sqrt(2*pi)*exp(erfinv(2*j-1)*(erfinv(2*j-1)+sqrt(2)))


for id=1:NumberOfElements MaterialParam[id]=exp(sqrt(2)*erfinv.(2 .* j .-1))*a end
solverparam=(elemtype =ElemType, Qpt=QuadPoints, Nelem=NumberOfElements, Order=parse(Int,ElemType[end]))
u1=solver1D.main(Nodes1,Elements,MaterialParam,solverparam)
u1=u1[Int64(floor(length(u1)/2))].*(exp(-(log(y))^2/2)/(y*sqrt(2*pi)))
output[ct]=u1*(1-0.00001-0+0.00001)
#output[ct]=(exp(-(log(y))^2/2)/(y*sqrt(2*pi)))
ct=ct+1
end

println(mean(output))
figure()
plot(eval_points,output,"*")

end


main()
