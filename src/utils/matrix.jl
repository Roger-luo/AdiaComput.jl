typealias BlasFloat Base.LinAlg.BlasFloat

#trotter expansion
function trotter(A::AbstractMatrix,B::AbstractMatrix,P::Int64)
    return (expm(full(A/(2*P)))*expm(full(B/P))*expm(full(A/(2*P))))^P
end

function (⊗)(A::AbstractMatrix,B::AbstractMatrix)
    @assert size(A)[1]==size(A)[2]
    @assert size(B)[1]==size(B)[2]

    return kron(A,B)
end

function (⊕)(A::AbstractSparseMatrix,B::AbstractSparseMatrix)
    return blkdiag(A,B)
end

function InsertDiagZero!(A::AbstractSparseMatrix)
    for i = 2:min(A.m,A.n)+1
        if A.colptr[i]-A.colptr[i-1]==0
            insert!(A.nzval,i-1,0)
            insert!(A.rowval,i-1,i-1)
            for j = i:min(A.m,A.n)+1
                A.colptr[j] += 1
            end
        end
    end
end

function diagexp(H::AbstractMatrix)
    return convert(Array{Complex128,1},pmap(exp,diag(H)))|>spdiagm
end
