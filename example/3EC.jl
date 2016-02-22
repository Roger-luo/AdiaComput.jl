using AdiaRoll
using PyPlot

blas_set_num_threads(4)
# eigen,state = evolution(1e3,[TruthTable(0b10,[1])],1)

eigen,state,p = evolution(
    1e4,
    [TruthTable(0b10010111,[1,2,3]),
    TruthTable(0b10010111,[2,3,4]),
    TruthTable(0b10010111,[3,4,5]),
    TruthTable(0b10010111,[5,6,7]),
    TruthTable(0b10010111,[1,5,6]),
    TruthTable(0b10010111,[2,3,6]), 
    ],7,dt=1)

# @show real(state)
State2TruthTable(state,3)
@show p



figure(1)
for i=1:size(eigen)[1]
    plot(real(eigen[i,:].'))
end

savefig("imag/3EC-7bit.png")
show()