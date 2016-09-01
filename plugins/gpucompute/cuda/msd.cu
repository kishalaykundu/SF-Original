#ifndef MSD_KERNEL_H
#define MSD_KERNEL_H

#include <cuda.h>
#include <cuda_runtime_api.h>

__constant__ float factor0;
__constant__ float factor1;
__constant__ unsigned int N;

texture<float4, cudaTextureType2D, cudaReadModeElementType> forceTexture;

// function to do time corrected Verlet integration for first two time-steps
extern "C" __global__ void displace_01 (float4 *src, float4 *dest, float2 *texCoords, float *mass)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;

  while (tid < N){

    float2 coords = texCoords [tid];
  	float4 force = tex2D (forceTexture, coords.x, coords.y);
  	float lmass_rec = mass [tid];

  	float4 acceleration = {lmass_rec*force.x, lmass_rec*force.y, lmass_rec*force.z, 1.};
  	float4 velocity = {factor0*acceleration.x, factor0*acceleration.y, factor0*acceleration.z, 1.};

    float4 prev = src [tid];
    float4 curr = {
                  prev.x + factor0*velocity.x + 0.5*factor1*acceleration.x,
                  prev.y + factor0*velocity.y + 0.5*factor1*acceleration.y,
                  prev.z + factor0*velocity.z + 0.5*factor1*acceleration.z,
                  1.
                  };
    dest [tid] = curr;

    tid += blockDim.x*gridDim.x;
  }
}

// function to do time corrected Verlet integration for time-steps 3 and later
extern "C" __global__ void displace_n (float4 *src, float4 *dest, float2 *texCoords, float *mass)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;

  while (tid < N){

    float2 coords = texCoords [tid];
  	float4 force = tex2D (forceTexture, coords.x, coords.y);
  	float lmass_rec = mass [tid];

  	float4 acceleration = {lmass_rec*force.x, lmass_rec*force.y, lmass_rec*force.z, 1.};

    float4 prev = dest [tid];
    float4 curr = src [tid];

    float4 future = {
                    curr.x + factor0*(curr.x - prev.x) + factor1*acceleration.x,
                    curr.y + factor0*(curr.y - prev.y) + factor1*acceleration.y,
                    curr.z + factor0*(curr.z - prev.z) + factor1*acceleration.z,
                    1.
                    };

    dest [tid] = future;

    tid += blockDim.x*gridDim.x;
  }
}

#endif
