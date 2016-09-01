# Top level file for the CU_MSD (CUDA accelerated eXtended Mass Spring Damper) module
project (SF_CUMSD CXX)

if (NOT WITH_VECTOR3 OR WITH_VECTOR3 STREQUAL "OFF")
    set (CUMSD_SRCS ${SF_SOURCE_DIR}/common/vec4.cpp)
endif ()

set (CUMSD_SRCS ${CUMSD_SRCS}
    ${SF_SOURCE_DIR}/plugins/gpucompute/cuda/msd.cu
    ${SF_SOURCE_DIR}/common/vec2.cpp
    ${SF_SOURCE_DIR}/common/vec3.cpp
    ${SF_SOURCE_DIR}/common/mat3x3.cpp
    ${SF_SOURCE_DIR}/common/GL/common.cpp
    ${SF_SOURCE_DIR}/common/GL/texture.cpp
    src/Common.cpp
    src/Mesh.cpp
    src/Plugin.cpp)

include_directories (./inc
    ${SF_PLUGINS_DIR}/physics/
    ${SF_SOURCE_DIR}/common/
    ${SF_ENGINE_DIR}/inc/
    ${SF_GRAPHICS_DIR}/Display/inc/
    /usr/include/libxml2/
    /usr/local/include/eigen2/)

# Find CUDA packages to be used
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/packages/FindCUDA/CMake/cuda" ${CMAKE_MODULE_PATH})
find_package (CUDA REQUIRED)

# Set compiler options for nvcc
set (${CUDA_NVCC_FLAGS} "-O3;-Wall")
cuda_wrap_srcs ("CudaMSD" PTX generated_ptx_files ${SF_SOURCE_DIR}/plugins/gpucompute/cuda/msd.cu)

set (CUMSD_LIBS ${CUMSD_LIBS} ${MATH_LIB} ${TIME_LIB} ${XML_LIB} ${BOOST_THREAD_LIB} ${NATIVE_THREAD_LIB} ${OPENGL_LIBRARY})

# Set name of the library
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (CUMSD_LIB "CudaMsd-debug")
else ()
    set (CUMSD_LIB "CudaMsd")
endif ()

# Final set up to build CUDA-enabled XFEM library
cuda_add_library (${CUMSD_LIB} SHARED ${CUMSD_SRCS})
target_link_libraries (${CUMSD_LIB} ${CUMSD_LIBS})

if (NOT CMAKE_BUILD_TYPE)
  set_target_properties (${CUMSD_LIB} PROPERTIES COMPILE_FLAGS "-std=c++0x -O2 -Wall")
else ()
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set_target_properties (${CUMSD_LIB} PROPERTIES COMPILE_FLAGS "-std=c++0x -DNDEBUG -O4 -Wall")
  else ()
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set_target_properties (${CUMSD_LIB} PROPERTIES COMPILE_FLAGS "-std=c++0x -g -fno-inline -Wall")
    endif ()
  endif ()
endif ()

# cuda_build_clean_target ()
