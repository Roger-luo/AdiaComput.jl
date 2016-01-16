#single problem Hamiltonian
function pcHamilton(truth_table::Integer,bitID::AbstractVector{Integer},bitnum::Integer)
    #bounds check
    @assert length(bitID)<=bitnum
    # @assert bitnum<length(bin(truth_table))

    temp = truth_table
    if length(bitID)<=bitnum
        temp=TruthTableExpand(temp,bitID,bitnum)
    end
    res = sparse(diagm([parse(Int,i) for i in temp]))
    return res
end

#final problem Hamiltonian
function pHamilton(expr::LogicExpr,bitnum::Integer)
    sum = pcHamilton(expr[1].value,expr[1].bitID,bitnum)

    for i = 2:length(expr)
        sum+=pcHamilton(expr[i].value,expr[1].bitID,bitnum)
    end
    return sum
end

function bHamilton(bitnum::Int)
    H = sparse(0.5*(I-σ₁))
    Iden = diagm([1,1])
    for i=2:bitnum
        H = H⊗Iden
    end

    for i=2:bitnum
        H_BC = Iden
        for j=2:bitnum
            if i==j
                H_BC = H_BC⊗sparse(0.5*(I-σ₁))
            else
                H_BC = H_BC⊗Iden
            end
        end
        H+=H_BC
    end
    return H
end

function Hamiltonian(
    s::Real,
    bHamiltonian::AbstractMatrix,
    pHamiltonian::AbstractMatrix
    )
    @assert 0<=s<=1
    return (1-s)*bHamiltonian+s*pHamiltonian
end