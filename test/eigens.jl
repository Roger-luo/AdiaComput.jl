using QuDynamics
using AdiaComput
using QuSAT
using QuBase
using QuComputStates

n = 5

ins,ans = generate(n)
pH = QuArray(pHamiltonian(ins,n),(comput_basis(n),comput_basis(n)))

aqc = AQC(pH,n;maxtime=1)
bH = bHamiltonian(aqc)

eigens = ((bH|>coeffs|>full|>eig)[1]|>real|>sort)[1:2]

for t = 0:1e-2:1
    H = Hamiltonian(aqc,t)
    append!(eigens,((H|>coeffs|>full|>eig)[1]|>real|>sort)[1:2])
end

eigens = reshape(eigens,2,100)

using PyPlot

figure()

eigens[1,:]|>plot
eigens[2,:]|>plot
