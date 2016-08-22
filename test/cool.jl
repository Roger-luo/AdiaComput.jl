using AdiaComput
using QuComputStates

s = BellState(2)

@show abs2(s[1])

H = [1 0 0 0;0 2 0 0;0 0 3 0;0 0 0 4]

gamma = 0
t     = 0.2

@assert -π/2 <= t-gamma <= π/2
@assert -π/2 <= 4t-gamma <= π/2

U = expm(-im*H*t)

@show size(U)

s = cool!(s,U,t,gamma)

@show abs2(s[1])
@show abs2(s[2])
@show abs2(s[3])
@show abs2(s[4])
