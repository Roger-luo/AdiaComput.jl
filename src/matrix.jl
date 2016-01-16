const σ₀=I
const σ₁=[0 1;1 0]
const σ₂=[0 -im;im 0]
const σ₃=[1 0;0 -1]


function (⊗)(A::AbstractMatrix,B::AbstractMatrix)
    @assert size(A)[1]==size(A)[2]
    @assert size(B)[1]==size(B)[2]

    return kron(A,B)
end

function (⊕)(A::AbstractMatrix,B::AbstractMatrix)
    return full(blkdiag(sparse(A),sparse(B)))
end