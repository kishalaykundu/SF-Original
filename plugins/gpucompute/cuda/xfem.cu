#ifndef XFEM_KERNEL_H
#define XFEM_KERNEL_H

#include <cuda.h>
#include <cuda_runtime_api.h>

extern "C" __global__ void conjugate_gradient (float *src, float *dest, unsigned int N)
{
  unsigned int ind = blockIdx.x*blockDim.x + threadIdx.x;
  if (ind < N){
    dest [ind] = src[ind] + 1.;
  }
}

#endif
