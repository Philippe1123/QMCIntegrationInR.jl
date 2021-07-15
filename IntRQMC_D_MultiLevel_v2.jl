

using Random, Distributions
using DigitalNets
using LatticeRules
using MultilevelEstimators
using PyPlot






function main()
    GC.gc()

#QMCType="Lattice"
QMCType="Sobol2"




    #Random.seed!(12345) # Setting the seed
    #d = MultilevelEstimators.Normal()
    #d = MultilevelEstimators.TruncatedNormal(0.0,1.0,-2.0,2.0)
    funTransform=(x,y)->transform(MultilevelEstimators.TruncatedNormal(0,1,-y,y),x)
#    funTransform=(x,y)->transform(MultilevelEstimators.Normal(),x)

f(x,gamma)=prod(1 .+ x .* gamma,dims=1)
distrb=Distributions.Normal()

b=-1:-0.5:-5.5
RequestedTolerance_vec=10 .^ b
RequestedTolerance_vec=vec(RequestedTolerance_vec)
time_vec=zeros(length(RequestedTolerance_vec))
err_vec=zeros(length(RequestedTolerance_vec))
truncationError_vec=zeros(length(RequestedTolerance_vec))
cubatureError_vec=zeros(length(RequestedTolerance_vec))

s=4
M=4
N=2  #2^N start number of samples
idx=1
while idx<=length(RequestedTolerance_vec)
GC.gc()


Matrx_std=zeros(length(N),1)
ct=1

level=12
sol=zeros(level+1)
sol_uncor=zeros(level+1)
box_fine=0
box_coarse=0
Sample_Multiplier=0
    l=0
    truncationError=10000
    cubatureError=10000
    soll=10000
RequestedTolerance=RequestedTolerance_vec[idx]
t=@elapsed begin
    percentagePtInSmallBox=[]
while l<=level



                p_f=2^(N+l+Sample_Multiplier)
                p_c=2^(N+l-1+Sample_Multiplier)


                p_f_FOR_BOX=2^(N+l)
                p_c_FOR_BOX=2^(N+l-1)

                println("---------")
                println("level ",l)
                println("Samples fine ",p_f)
                println("Samples coarse ",p_c)


                Matrix_f=zeros(s,p_f,M)
                Matrix_c=zeros(s,p_f,M)
                sizeShift=zeros(M)
                a=1:1:s
                gamma=1 ./ a .^2

                if(QMCType=="Lattice")
                    Lattice=LatticeRule(s)

                elseif(QMCType=="Sobol2")
                    Lattice=DigitalNet64(s)
                end
            #
                    for Nshift=1:M

                     if(QMCType=="Lattice")
                         shiftLat=ShiftedLatticeRule(Lattice)

                     elseif(QMCType=="Sobol2")
                         shiftLat=DigitalShiftedDigitalNets64(Lattice)

                     end




                        ct=1
                            for id=1:length(shiftLat[0:p_f_FOR_BOX])-1
                                box_fine=sqrt(2*2*log(p_f_FOR_BOX))
                                box_coarse=sqrt(2*2*log(p_c_FOR_BOX))

                                Matrix_f[:,id,Nshift]=map.(funTransform,shiftLat[id-1],box_fine)

                                cumsum=sum((Matrix_f[:,id,Nshift] .> box_coarse)  .| (Matrix_f[:,id,Nshift] .< -box_coarse))
                                if(cumsum>0)
                                else
                                Matrix_c[:,ct,Nshift]=Matrix_f[:,id,Nshift]
                                ct=ct+1
                                end
                            end
                        sizeShift[Nshift]=ct-1

                    end
            #       println(sizeShift)
                  # println(Matrix_c)

                    Matrix_c=Matrix_c[:,1:Int64(minimum(sizeShift)),:]
                    nshifts_effective=minimum(sizeShift)
                    pct=(1-size(Matrix_c,2)/size(Matrix_f,2))*100
                    println("Percentage of points in smaller box is ",pct , " %")
                    println("box fine is ",-box_fine," ",box_fine)
                    println("box coarse is ",-box_coarse," ",box_coarse)

                    append!(percentagePtInSmallBox,pct)

                G_fine=mean(f(Matrix_f,gamma),dims=2)
            #    println(G_fine)
                corrfactor_fine=cdf(distrb,box_fine)-cdf(distrb,-box_fine)
                corrfactor_fine.^s
                G_coarse=mean(f(Matrix_c,gamma),dims=2)
                corrfactor_coarse=cdf(distrb,box_coarse)-cdf(distrb,-box_coarse)
                corrfactor_coarse=corrfactor_coarse.^s
            #    println(G_coarse)
                diff=(G_fine .- G_coarse)*(corrfactor_fine-corrfactor_coarse)
                println("corrfactor big box ",corrfactor_fine)
                println("corrfactor small box ",corrfactor_coarse)
                println("corrfactor diff ",(corrfactor_fine-corrfactor_coarse))



                #QMc_R=mean(diff,dims=2)
                QMc_Q=mean((diff),dims=3)
                sol[l+1]=QMc_Q[1]

                #sol_uncor[l+1]=mean((G_fine .- G_coarse),dims=3)[1]



                QMc_R=mean(G_fine.*corrfactor_fine,dims=2)
                QMc_std=std(QMc_R)/sqrt(nshifts_effective)


                cubatureError=QMc_std
                truncationError=sol[l+1]


                println("Cubature error is ",cubatureError)
                println("Truncation error is ",sol[l+1])




                soll=mean((G_fine*corrfactor_fine),dims=3)[1]
                println("Sol is ",soll)
                println("abs error is ", abs(1-soll))





                if(cubatureError>RequestedTolerance*0.5 || truncationError==0)
                    Sample_Multiplier=Sample_Multiplier+1
                else
                    if(abs(truncationError)>RequestedTolerance*0.5 && truncationError!=0)
                        l=l+1
                    else
                        l=level+2
                       truncationError_vec[idx]=truncationError
                       cubatureError_vec[idx]=cubatureError



                    end
                end




end

figure()
plot(percentagePtInSmallBox,"-*")

end
println("Runtime is ",t," sec")
println("******************************************************************************")

time_vec[idx]=t
err_vec[idx]=abs(1-soll)
idx=idx+1
end

figure()
loglog(RequestedTolerance_vec,time_vec,"-*")
title("Runtime in function of requested tolerance")
xlabel("RMSE")
ylabel("run time sec")
figure()
loglog(RequestedTolerance_vec,err_vec,"-*")
println(time_vec)
loglog(RequestedTolerance_vec,cubatureError_vec,"-*")
loglog(RequestedTolerance_vec,abs.(truncationError_vec),"-*")

str=string(QMCType, "\n" ,"s = " , s, "shift = ",M ," \n","Errors in function of requested RMSE")

title(str)
ylabel("Abs. error")
xlabel("RMSE")
legend(("abs. error on sol","cubature error","abs. trunction error"))
println(truncationError_vec)
println(cubatureError_vec)
GC.gc()
#println(sol)
#println(sum(sol))
#println(sol_uncor)
#println(sum(sol_uncor))
end


main()
