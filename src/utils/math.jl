include("matrix.jl")

import Base.LinAlg.norm

function norm(vec::CudaArray)
    return CUBLAS.nrm2(vec)
end

function normalize!(vec::CudaArray)
    CUBLAS.scal!(vec.dims[1],1/norm(vec),vec,1)
end

function normalize!(vec::AbstractVector)
    return vec[:]=vec./norm(vec)
end

export normalize!
