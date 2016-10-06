module AdiaComput
using QuComputStates,QuDynamics,JuMP,QuSAT,QuBase,Ipopt, ExpmV, Expokit

import Base: |>
import QuDynamics: QuStateEvolution,operator,propagate
import QuBase: AbstractQuMatrix,AbstractQuVector,similar_type
export AQC,bHamiltonian,pHamiltonian,Hamiltonian,|>,operator

include("Hamiltonian.jl")
include("AQCShrodingerEq.jl")
include("AlgorithmCooling.jl")

typealias AbstractQuVecOrMat Union{AbstractQuVector,AbstractQuMatrix}

type AQC{N,H<:AbstractQuMatrix}<:AbstractQC{N}
    eq::AQCShrodingerEq
    state::AbstractQuVecOrMat
    tlist
    method

    function AQC(HP::H,tlist,maxtime,method)
        @assert tlist[end] <= maxtime "invalid end time of time list"

        eq = AQCShrodingerEq(HP,maxtime)
        # tlist = 0.:dt:maxtime
        state = QuArray(convert(Array{Complex128,1},[1/sqrt(2^N) for i=1:2^N]),comput_basis(N))
        new(eq,state,tlist,method)
    end
end

# basic simulation on default timeline
AQC{H<:AbstractQuMatrix}(HP::H,n::Int;method = QuODE45(),dt = 1e-2, maxtime = 1.0) = AQC{n,H}(HP,0.:dt:maxtime,maxtime,method)

# modified timeline
AQC{H<:AbstractQuMatrix}(HP::H,n::Int,tlist;method = QuODE45(), maxtime = 1.0) = AQC{n,H}(HP,tlist,maxtime, method)


# Hamiltonian accessor
bHamiltonian(aqc::AQC) = aqc.eq.HB
pHamiltonian(aqc::AQC) = aqc.eq.HP

Hamiltonian(aqc::AQC,t) = operator(aqc.eq,t)

function (|>)(tlist::Range,aqc::AQC)
    aqc.tlist = tlist
end

function (|>)(procedure,aqc::AQC)
    aqc.tlist = procedure[1]
    aqc.method = procedure[2]
end

function Base.start(aqc::AQC)
    init_state = aqc.state
    t_state = start(aqc.tlist)
    t,t_state = next(aqc.tlist,t_state)
    return QuPropagatorState(init_state,t,t_state)
end

function Base.next(aqc::AQC, qustate::QuPropagatorState)
    current_qustate = qustate.state
    current_t = qustate.t
    t,t_state = next(aqc.tlist, qustate.t_state)
    next_qustate = aqc.state = propagate(aqc.method, aqc.eq, t, current_t, current_qustate)
    return (t, next_qustate), QuPropagatorState(next_qustate, t, t_state)
end

Base.done(aqc::AQC, qustate::QuPropagatorState) = done(aqc.tlist, qustate.t_state)
end
