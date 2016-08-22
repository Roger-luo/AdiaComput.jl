using QuDynamics
using AdiaComput
using QuSAT
using QuBase
using QuComputStates


function aqc_test(n;dt=1e-2)
    p = 0
    for k = 1:100
        ins,ans = generate(n)
        pH = QuArray(pHamiltonian(ins,n),(comput_basis(n),comput_basis(n)))

        aqc = AQC(pH,n;maxtime=1,method=QuExpokit())
        bH = bHamiltonian(aqc)

        # START

        (0.:1:10,DemonCooling(coeffs(pH)))|>aqc
        for (t,psi) in aqc end

        p += norm(aqc.state[ans[1]+1])^2
    end
    return p/100
end

@show aqc_test(4)

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


# temp = (coeffs(0.4*bH+0.6*pH)|>full|>eig)
# @show (temp[1]|>real|>findmin)[1]
# @show (temp[1]|>real|>findmax)[1]
# id = (temp[1]|>real|>findmin)[2]
# ground = temp[2][:,id]

# @show norm(coeffs(aqc.state).'*ground)
