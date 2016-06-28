# using CUDArt,CUSPARSE

import Base:+,-,*,/,.*,.+,.-,./
export diagexp

rt = CUDArt.CUDArt_gen

function diagexp(A::CudaArray{Complex64})
    md = CuModule("$dir/src/utils/cuda/cuMatrix.ptx", false)
    diagexp = CuFunction(md, "diagexp_cf")
    nsm = attribute(device(), rt.cudaDevAttrMultiProcessorCount)
    mul = min(32, ceil(Int, length(A)/(256*nsm)))
    expH = CudaArray(Complex64,size(A)...)
    launch(diagexp, mul*nsm, 256, (A,expH,length(A)))
    return expH
end

function diagexp(A::CudaArray{Complex128};devlist=[0])
    md = CuModule("$dir/src/utils/cuda/cuMatrix.ptx", false)
    diagexp = CuFunction(md, "diagexp_df")
    nsm = attribute(device(), rt.cudaDevAttrMultiProcessorCount)
    mul = min(32, ceil(Int, length(A)/(256*nsm)))
    expH = CudaArray(Complex128,size(A)...)
    launch(diagexp, mul*nsm, 256, (A,expH,length(A)))
    return expH
end

function diagexp{T}(A::CudaSparseMatrixCSR{T})
    A.nzVal = diagexp(A.nzVal)
    return A
end

function diagexp{T}(A::CudaSparseMatrixCSC{T})
    A.nzVal = diagexp(A.nzVal)
    return A
end

##############################
#  CUBLAS operator wrappers
##############################

function (+){T}(a::CudaArray{T},b::CudaArray{T})
    c = CudaArray(T,b.dims[1])
    CUBLAS.copy!(c,b)
    CUBLAS.axpy!(one(T),a,c)
    return c
end

function (.+){T}(a::CudaArray{T},b::CudaArray{T})
    CUBLAS.axpy!(one(T),b,a)
    return a
end

function (-){T}(a::CudaArray{T},b::CudaArray{T})
    c = CudaArray(T,a.dims[1])
    CUBLAS.copy!(c,a)
    CUBLAS.axpy!(-one(T),b,c)
    return c
end

function (.-){T}(a::CudaArray{T},b::CudaArray{T})
    CUBLAS.axpy!(-one(T),b,a)
    return a
end

function (*){T}(alpha::T,a::CudaArray{T})
    c = CudaArray(T,a.dims[1])
    CUBLAS.copy!(c,a)
    CUBLAS.scal!(c.dims[1],alpha,c,1)
    return c
end

(*){T}(a::CudaArray{T},alpha::T) = (*)(alpha,a)

function (.*){T}(a::CudaArray{T},alpha::T)
    CUBLAS.scal!(a.dims[1],alpha,a,1)
    return a
end

(.*){T}(alpha::T,a::CudaArray{T}) = (.*)(a,alpha)

(/){T}(a::CudaArray{T},alpha::T) = (*)(1/alpha,a)
(/){T}(alpha::T,a::CudaArray{T}) = (*)(1/alpha,a)

(./){T}(a::CudaArray{T},alpha::T) = (.*){T}(a,1/alpha)
(./){T}(alpha::T,a::CudaArray{T}) = (.*){T}(a,1/alpha)



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

function diagexp2(H::CudaSparseMatrix)
    return convert(Array{Complex128,1},pmap(exp,diag(to_host(H))))|>spdiagm|>CudaSparseMatrixCSR
end
