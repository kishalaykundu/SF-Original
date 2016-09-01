/**
 * @file common.h
 * @author Kishalay Kundu <kishalay.kundu@gmail.com>
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 * Common functions used by CUDA based functions.
 */
#pragma once

#include <cassert>
#include <cstdio>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

extern "C"{
#include <cuda.h>
#include <cuda_runtime_api.h>
}

namespace SF {

  // function to check for any CUDA error
  inline void checkCUDAErrorPrivate (cudaError_t &error, const string &file, unsigned int line)
  {
    if ((error = cudaGetLastError ()) != cudaSuccess){
        fprintf ( stdout, "%s[%u]:\tCUDA error: %s\n", basename (file.c_str ()), line, cudaGetErrorString (error));
    }
  }
  #define checkCUDAError(error) (checkCUDAErrorPrivate (error, __FILE__, __LINE__))

  // function to check for CUresult
  void printCUResultPrivate (CUresult &result, const string &file, unsigned int line)
  {
#ifndef NDEBUG
    if (result == CUDA_SUCCESS){
      return;
    }
    char code [64];
    switch (result){
      case CUDA_ERROR_INVALID_VALUE:
        sprintf (code, "CUDA_ERROR_INVALID_VALUE");
        break;
      case CUDA_ERROR_OUT_OF_MEMORY:
        sprintf (code, "CUDA_ERROR_OUT_OF_MEMORY");
        break;
      case CUDA_ERROR_NOT_INITIALIZED:
        sprintf (code, "CUDA_ERROR_NOT_INITIALIZED");
        break;
      case CUDA_ERROR_DEINITIALIZED:
        sprintf (code, "CUDA_ERROR_DEINITIALIZED");
        break;
      case CUDA_ERROR_PROFILER_DISABLED:
        sprintf (code, "CUDA_ERROR_PROFILER_DISABLED");
        break;
      case CUDA_ERROR_PROFILER_NOT_INITIALIZED:
        sprintf (code, "CUDA_ERROR_PROFILER_NOT_INITIALIZED");
        break;
      case CUDA_ERROR_PROFILER_ALREADY_STARTED:
        sprintf (code, "CUDA_ERROR_PROFILER_ALREADY_STARTED");
        break;
      case CUDA_ERROR_PROFILER_ALREADY_STOPPED:
        sprintf (code, "CUDA_ERROR_PROFILER_ALREADY_STOPPED");
        break;
      case CUDA_ERROR_NO_DEVICE:
        sprintf (code, "CUDA_ERROR_NO_DEVICE");
        break;
      case CUDA_ERROR_INVALID_DEVICE:
        sprintf (code, "CUDA_ERROR_INVALID_DEVICE");
        break;
      case CUDA_ERROR_INVALID_IMAGE:
        sprintf (code, "CUDA_ERROR_INVALID_IMAGE");
        break;
      case CUDA_ERROR_INVALID_CONTEXT:
        sprintf (code, "CUDA_ERROR_INVALID_CONTEXT");
        break;
      case CUDA_ERROR_CONTEXT_ALREADY_CURRENT:
        sprintf (code, "CUDA_ERROR_CONTEXT_ALREADY_CURRENT");
        break;
      case CUDA_ERROR_MAP_FAILED:
        sprintf (code, "CUDA_ERROR_MAP_FAILED");
        break;
      case CUDA_ERROR_UNMAP_FAILED:
        sprintf (code, "CUDA_ERROR_UNMAP_FAILED");
        break;
      case CUDA_ERROR_ARRAY_IS_MAPPED:
        sprintf (code, "CUDA_ERROR_ARRAY_IS_MAPPED");
        break;
      case CUDA_ERROR_ALREADY_MAPPED:
        sprintf (code, "CUDA_ERROR_ALREADY_MAPPED");
        break;
      case CUDA_ERROR_NO_BINARY_FOR_GPU:
        sprintf (code, "CUDA_ERROR_NO_BINARY_FOR_GPU");
        break;
      case CUDA_ERROR_ALREADY_ACQUIRED:
        sprintf (code, "CUDA_ERROR_ALREADY_ACQUIRED");
        break;
      case CUDA_ERROR_NOT_MAPPED:
        sprintf (code, "CUDA_ERROR_NOT_MAPPED");
        break;
      case CUDA_ERROR_NOT_MAPPED_AS_ARRAY:
        sprintf (code, "CUDA_ERROR_NOT_MAPPED_AS_ARRAY");
        break;
      case CUDA_ERROR_NOT_MAPPED_AS_POINTER:
        sprintf (code, "CUDA_ERROR_NOT_MAPPED_AS_POINTER");
        break;
      case CUDA_ERROR_ECC_UNCORRECTABLE:
        sprintf (code, "CUDA_ERROR_ECC_UNCORRECTABLE");
        break;
      case CUDA_ERROR_UNSUPPORTED_LIMIT:
        sprintf (code, "CUDA_ERROR_UNSUPPORTED_LIMIT");
        break;
      case CUDA_ERROR_CONTEXT_ALREADY_IN_USE:
        sprintf (code, "CUDA_ERROR_CONTEXT_ALREADY_IN_USE");
        break;
      case CUDA_ERROR_INVALID_SOURCE:
        sprintf (code, "CUDA_ERROR_INVALID_SOURCE");
        break;
      case CUDA_ERROR_FILE_NOT_FOUND:
        sprintf (code, "CUDA_ERROR_FILE_NOT_FOUND");
        break;
      case CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND:
        sprintf (code, "CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND");
        break;
      case CUDA_ERROR_SHARED_OBJECT_INIT_FAILED:
        sprintf (code, "CUDA_ERROR_SHARED_OBJECT_INIT_FAILED");
        break;
      case CUDA_ERROR_OPERATING_SYSTEM:
        sprintf (code, "CUDA_ERROR_OPERATING_SYSTEM");
        break;
      case CUDA_ERROR_INVALID_HANDLE:
        sprintf (code, "CUDA_ERROR_INVALID_HANDLE");
        break;
      case CUDA_ERROR_NOT_FOUND:
        sprintf (code, "CUDA_ERROR_NOT_FOUND");
        break;
      case CUDA_ERROR_NOT_READY:
        sprintf (code, "CUDA_ERROR_NOT_READY");
        break;
      case CUDA_ERROR_LAUNCH_FAILED:
        sprintf (code, "CUDA_ERROR_LAUNCH_FAILED");
        break;
      case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES:
        sprintf (code, "CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES");
        break;
      case CUDA_ERROR_LAUNCH_TIMEOUT:
        sprintf (code, "CUDA_ERROR_LAUNCH_TIMEOUT");
        break;
      case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING:
        sprintf (code, "CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING");
        break;
      case CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED:
        sprintf (code, "CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED");
        break;
      case CUDA_ERROR_PEER_ACCESS_NOT_ENABLED:
        sprintf (code, "CUDA_ERROR_PEER_ACCESS_NOT_ENABLED");
        break;
      case CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE:
        sprintf (code, "CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE");
        break;
      case CUDA_ERROR_CONTEXT_IS_DESTROYED:
        sprintf (code, "CUDA_ERROR_CONTEXT_IS_DESTROYED");
        break;
      case CUDA_ERROR_ASSERT:
        sprintf (code, "CUDA_ERROR_ASSERT");
        break;
      case CUDA_ERROR_TOO_MANY_PEERS:
        sprintf (code, "CUDA_ERROR_TOO_MANY_PEERS");
        break;
      case CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED:
        sprintf (code, "CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED");
        break;
      case CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED:
        sprintf (code, "CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED");
        break;
      case CUDA_ERROR_UNKNOWN:
        sprintf (code, "CUDA_ERROR_UNKNOWN");
        break;

      default:
        sprintf (code, "CUDA_UNKNOWN_RESULT");
    }
    fprintf (stdout, "%s[%u]: %s\n", basename (file.c_str ()), line, code);
#endif
  }
  #define printCUResult(status) (printCUResultPrivate (status, __FILE__, __LINE__))

  // function to read CUDA PTX file
  inline void readCudaPTXFile (const string &file, string &contents)
  {
    assert (!file.empty ());
    ifstream is;
    is.open (file.c_str (), ios::binary);
    assert (!is.fail ());
    is.seekg (0, ios::end);
    size_t length = is.tellg ();
    is.seekg (0, ios::beg);
    char *src = new char [length + 1];
    is.read (src, length);
    is.close ();
    src [length] = '\0';

    contents = string (src);
  }
}
