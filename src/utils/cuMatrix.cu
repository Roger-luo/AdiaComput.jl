#include <stdint.h>
#include <math.h>

typedef struct
{
    double x;
    double y;
}complexDouble;

typedef struct
{
    float x;
    float y;
}complex;

template <typename T>
__device__ void diagexp(T *A,T *C, int len) {
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if(idx<len)
    {
        C[idx].x = exp(A[idx].x)*cos(A[idx].y);
        C[idx].y = exp(A[idx].x)*sin(A[idx].y);
    }
}

extern "C"
{
    void __global__ diagexp_cf(complex *A,complex *C, int len){diagexp(A,C,len);};
    void __global__ diagexp_df(complexDouble *A,complexDouble *C, int len){diagexp(A,C,len);};
}
