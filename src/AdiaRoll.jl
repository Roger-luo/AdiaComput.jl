module AdiaRoll

import Base:bin

include("truthtable.jl")
include("matrix.jl")
include("Hamiltonian.jl")


#adiabatic evolution for one time step,a simple Euler form is used at present
function next_timestep(
    state::AbstractVector,
    curstep::Int64,
    bHamiltonian::AbstractMatrix,
    pHamiltonian::AbstractMatrix,
    evotime::Real,
    dt=1e-2
    )
    @show typeof(-im*dt*full(Hamiltonian(curstep*dt/evotime,bHamiltonian,pHamiltonian)))
    return expm(-im*dt*full(Hamiltonian(curstep*dt/evotime,bHamiltonian,pHamiltonian)))*state
end


doc"""
```julia
function evolution(evotime::Real,expr::LogicExpr,bitnum::Int;dt=1)->eigenvalue,state,success_probablity
```
`evolution` evolutes the adiabatic system

Parameters
---
- `evotime::Real` evolution time
- `expr::LogicExpr` LogicExpression,a collection of TruthTable
- `bitnum::Int` number of qubits
- `dt` time step,defaut is 1

Returns
---
- eigenvalue
- state
- success probablity defined as $|<\phi|\psi(T)>|^2$,where $\phi$ is the ground state of problem Hamiltonian,$\psi$ is the final state of the evolution
"""
function evolution(
    evotime::Real,
    expr::LogicExpr,
    bitnum::Int;
    dt=1
    )
    state = [1/sqrt(2^bitnum) for i=1:2^bitnum]
    bHamiltonian = bHamilton(bitnum)
    pHamiltonian = pHamilton(expr,bitnum)

    ##start calculation
    # eigenvalue = collect(complex(conj(state).'*Hamiltonian(0,bHamiltonian,pHamiltonian)*state))
    eigentemp   = eig(full(Hamiltonian(0,bHamiltonian,pHamiltonian)))
    eigenvalue  = eigentemp[1]
    eigennum    = length(eigenvalue)

    for i=1:dt:evotime/dt
        state = next_timestep(state,Int64(i),bHamiltonian,pHamiltonian,evotime,dt)
        # append!(eigenvalue, conj(state).'*Hamiltonian(i*dt/evotime,bHamiltonian,pHamiltonian)*state)
        eigentemp = eig(full(Hamiltonian(i/evotime,bHamiltonian,pHamiltonian)))
        append!(eigenvalue, eigentemp[1])
    end

    success_prob= abs2([x==findmin(diag(pHamiltonian))[2]?1:0 for x=1:2^bitnum].'*state)
    return reshape(eigenvalue,eigennum,Int(length(eigenvalue)/eigennum)),state,success_prob[1]
end

export evolution,LogicExpr,TruthTable,State2TruthTable
end
