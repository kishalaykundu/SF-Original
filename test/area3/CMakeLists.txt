# CMake file for test-area1
project (SF_TEST3 CXX)

# Set source file names
set (TEST3_SRCS
    ${SF_SOURCE_DIR}/plugins/gpucompute/cuda/xfem.cu
    ${SF_SOURCE_DIR}/common/vec4.cpp
    src/main.cpp)

include_directories (${SF_SOURCE_DIR}/common
                    ${SF_GRAPHICS_DIR}/Display/inc
                    ./inc)

# Find CUDA packages to be used
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/packages/FindCUDA/CMake/cuda" ${CMAKE_MODULE_PATH})
find_package (CUDA REQUIRED)

# Set compiler options for nvcc
set (${CUDA_NVCC_FLAGS} "-O3;-Wall")
cuda_wrap_srcs ("CudaXfem" PTX generated_ptx_files ${SF_SOURCE_DIR}/plugins/gpucompute/cuda/xfem.cu)

set (TEST3_LIBS ${BOOST_THREAD_LIB} ${MATH_LIB} ${NATIVE_THREAD_LIB} ${OPENGL_LIBRARY} ${GLU_LIBRARY} ${GLUT_LIB})

# Set final name of executable
cuda_add_executable (test3 ${TEST3_SRCS})
target_link_libraries (test3 ${TEST3_LIBS})

if (NOT CMAKE_BUILD_TYPE)
  set_target_properties (test3 PROPERTIES COMPILE_FLAGS "-std=c++0x -O2 -Wall")
else ()
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set_target_properties (test3 PROPERTIES COMPILE_FLAGS "-std=c++0x -DNDEBUG -O4 -Wall")
  else ()
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set_target_properties (test3 PROPERTIES COMPILE_FLAGS "-std=c++0x -g -fno-inline -Wall")
    endif ()
  endif ()
endif ()
