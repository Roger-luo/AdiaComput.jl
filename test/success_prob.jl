using QuDynamics
using AdiaComput
using QuSAT
using QuBase
using QuComputStates

function aqc_test(n,div;dt=1e-2)
    res = Float64[]

    ins,ans = generate(n)
    pH = QuArray(pHamiltonian(ins,n),(comput_basis(n),comput_basis(n)))

    aqc = AQC(pH,n;maxtime=1,method=QuExpokit())
    bH = bHamiltonian(aqc)

    # START
    (0.0:dt:div[1],QuODE45())|>aqc
    for (t,psi) in aqc end

    temp = (coeffs((1-div[1])*bH+div[1]*pH)|>full|>eig)
    @show (temp[1]|>real|>findmin)[1]
    @show (temp[1]|>real|>findmax)[1]
    id = (temp[1]|>real|>findmin)[2]
    ground = temp[2][:,id]

    push!(res,norm(coeffs(aqc.state).'*ground))


    for i = 1:length(div)-1
        (div[i]:dt:div[i+1],QuODE45())|>aqc
        for (t,psi) in aqc end

        temp = (coeffs((1-div[i+1])*bH+div[i+1]*pH)|>full|>eig)
        @show (temp[1]|>real|>findmin)[1]
        @show (temp[1]|>real|>findmax)[1]
        id = (temp[1]|>real|>findmin)[2]
        ground = temp[2][:,id]

        push!(res,norm(coeffs(aqc.state).'*ground))
    end

    (div[end]:dt:1,QuODE45())|>aqc
    for (t,psi) in aqc end

    temp = (coeffs(pH)|>full|>eig)
    @show (temp[1]|>real|>findmin)[1]
    @show (temp[1]|>real|>findmax)[1]
    id = (temp[1]|>real|>findmin)[2]
    ground = temp[2][:,id]

    push!(res,norm(coeffs(aqc.state).'*ground))


    return res
end

p = aqc_test(4,0.1:0.1:0.9)

using PyPlot

figure()
plot(p)
savefig("ground.png")

# res = Float64[]
# for i = 0.1:0.1:0.8
#     push!(res,aqc_test(4,(i,i+0.1)))
# end

# using PyPlot
#
# plot(res)
#
# savefig("res.png")

#
# 4 bits
# DATA
# normal
# p: 0.07244126481044265
# 0:0.6:0.8:1.0
# p: 0.6745965055249651
# 0:0.2:0.4:1.0
# p: 0.13587412610882968
# 0:0.8:0.9:1.0
# p: 0.9112006362656903
