module AdiaComput

using CUDArt,CUSPARSE,CUBLAS

abstract QuComput

const dir = Pkg.dir("AdiaComput")

include("utils.jl")
include("Base.jl")
include("Hamiltonian.jl")
include("Operator.jl")

end # module
