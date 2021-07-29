using Random, Distributions
using DigitalNets
using LatticeRules
using MultilevelEstimators
using PyPlot




function main()


id=0:1:20
nbPt= 2 .^ id
nbPt=1:10000:2^20
s=1

Lattice=ShiftedLatticeRule(LatticeRule(s))
DigitalNet=DigitalShiftedDigitalNets64(DigitalNet64(s))
Digital32=DigitalShiftedDigitalNets32(DigitalNet32(s))
time_lat=zeros(length(nbPt),1)
time_net=zeros(length(nbPt),1)
time_net32=zeros(length(nbPt),1)

ct=1
for i in nbPt

t_lat=@elapsed begin
    Lattice[0:i-1]
end
time_lat[ct]=t_lat

t_net=@elapsed begin
    DigitalNet[0:i-1]
end
time_net[ct]=t_net


t_net32=@elapsed begin
    Digital32[0:i-1]
end
time_net32[ct]=t_net32

ct=ct+1
println(i)
end

figure()
plot(nbPt,time_lat,"-*")
plot(nbPt,time_net,"-*")
plot(nbPt,time_net32,"-*")
legend(("Lattice","DigitalNet interlace 2","DigitalNet"))
xlabel("Nb. Points")
ylabel("Time [sec]")
title("Timings Shifted QMC")

figure()
plot(nbPt,time_net32./time_lat,"-*")
title("Timings Difference Shifted QMC")

xlabel("Nb. Points")
ylabel("Time [sec]")

end
main()
