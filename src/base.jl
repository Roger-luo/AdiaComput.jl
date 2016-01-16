abstract AbstractQuArray{T,N}
abstract AbstractOperator

typealias AbstractBaseExpr Union{AbstractVector{Function},Void}

type QuArray{T<:Number,N}<:AbstractQuArray{T,N}
    coeffs::AbstractArray{T,N}
    bases::AbstractVector
    baseExprs::AbstractBaseExpr

    function QuArray(coeffs::AbstractArray{T,N},bases::AbstractVector,baseExprs::AbstractBaseExpr)
        if typeof(baseExprs)==Void
            if N==1&&length(coeffs)==length(bases)
                new(coeffs,bases,baseExprs)
            elseif N==2&&size(coeffs)[1]==size(coeffs)[2]==length(bases)
                new(coeffs,bases,baseExprs)
            elseif typeof(baseExprs)==Void
                new(coeffs,bases)
            end
        elseif (N==1&&length(coeffs)==length(bases)==length(baseExprs))||(N==2&&size(coeffs)[1]==size(coeffs)[2]==length(bases)==length(baseExprs))
                new(coeffs,bases,baseExprs)
        else
            error("size do not match!\n")
        end
    end
end

QuArray{T<:Number,N}(coeffs::AbstractArray{T,N},bases::AbstractVector,baseExprs::AbstractBaseExpr)=QuArray{T,N}(coeffs,bases,baseExprs)
QuArray{T,N}(coeffs::AbstractArray{T,N},bases::AbstractVector)=QuArray(coeffs,bases,nothing)
QuArray{T,N}(coeffs::AbstractArray{T,N})=QuArray(coeffs,[1:size(coeffs)[1]],nothing)


function show(WaveFunc::QuArray)
    n = length(WaveFunc.coeffs)

    print("$(WaveFunc.coeffs[1])|$(WaveFunc.bases[1])⟩")
    for i = 2:n
        if WaveFunc.coeffs[i]>0
            print("+$(WaveFunc.coeffs[i])|$(WaveFunc.bases[i])⟩")
        elseif WaveFunc.coeffs[i]<0
            print("$(WaveFunc.coeffs[i])|$(WaveFunc.bases[i])⟩")
        end
        print("\n")
    end

    if typeof(WaveFunc.baseExprs)==Void
        print("[$(typeof(WaveFunc.coeffs)),$(size(WaveFunc.coeffs))]⊗[AbstractBase,$(size(WaveFunc.coeffs)[1])]\n")
    else
        print("[$(typeof(WaveFunc.coeffs)),$(size(WaveFunc.coeffs))]⊗[ExprBase,$(length(WaveFunc.baseExprs))]\n")
    end
end


macro show(WaveFunc)
    return :(show($WaveFunc))
end

# base(x,n)=x->sin(n*x)

# test = QuArray([1,1],[-1//2,1//2],[x->base(x,1)])
# @show test