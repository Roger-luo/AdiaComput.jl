function eigenvector(index::Integer,n::Integer)
    res = 1
    for i in bin(index,n)
        if i=='0'
            res = kron(res,1/sqrt(2)*[1,1])
        elseif i=='1'
            res = kron(res,1/sqrt(2)*[1,-1])
        end
    end
    return res
end

type AdiaComputer <: QuComput
    ###########################################
    #  Eternal
    ###########################################
    HB::AbstractMatrix      # Base Hamiltonian
    HP::AbstractMatrix      # Problem Hamiltonian
    P::AbstractMatrix       # SimilarMatrix for HB
    invP::AbstractMatrix    # inv(P)
    maxtime::Real           # max evolution time
    n::Int64                # number of bits
    dt::Real                # time step

    ###########################################
    # variable for current state
    ###########################################
    location::Real          # time location
    state::Union{AbstractVector,AbstractCudaArray}
    # eigens::AbstractVecOrMat

    ###########################################
    # measure
    ###########################################

    prob::Real

    ###########################################
    # options
    ###########################################
    GPU::Bool

    function AdiaComputer{M,N}(
        ins::Instance{M,N},
        n::Int, # number of bits
        maxtime::Real # max evolution time
        ;
        dt=1e-2,
        GPU::Bool = false
        )

        HB = bHamilton(n)
        HP = pHamilton(ins,n)
        # Similar Matrix
        P = eigenvector(0,n)
        for i=1:2^n-1
            P = [P eigenvector(i,n)]
        end
        invP = inv(P)
        # prepare the initial state
        state = convert(Array{Complex128,1},[1/sqrt(2^n) for i=1:2^n])
        # probility
        prob = 1
        # set time location to be the beginning
        location = 0.0

        # options

        if GPU
            InsertDiagZero!(HB)
            InsertDiagZero!(HP)

            HB = CudaSparseMatrixCSR(HB)
            HP = CudaSparseMatrixCSR(HP)
            global expH = CudaArray(zeros(Complex128,2^n))
            P = convert(Array{Complex128,2},P)
            invP = convert(Array{Complex128,2},invP)
            P  = CudaSparseMatrixCSR(sparse(P))
            invP = CudaSparseMatrixCSR(sparse(invP))
            state = CudaArray(state)
        end
        new(HB,HP,P,invP,maxtime,n,dt,location,state,prob,GPU)
    end
end

function Hamiltonian(Hs::AdiaComputer)
    return Complex128(1.0-Hs.location)*Hs.HB+Complex128(Hs.location)*Hs.HP
end

function reset!(Hs::AdiaComputer)
    Hs.state = convert(Array{Complex128,1},[1/sqrt(2^Hs.n) for i=1:2^Hs.n])
    if Hs.GPU
        Hs.state = CudaArray(Hs.state)
    end
    Hs.location = 0.0
    Hs.prob = 1
end

export AdiaComputer,reset!
