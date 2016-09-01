# Top level file for the CU_MSD (CUDA accelerated eXtended Mass Spring Damper) module
project (SF_CPUMSD CXX)

if (NOT WITH_VECTOR3 OR WITH_VECTOR3 STREQUAL "OFF")
    set (CPUMSD_SRCS ${SF_SOURCE_DIR}/common/vec4.cpp)
endif ()

set (CPUMSD_SRCS ${CPUMSD_SRCS}
    ${SF_SOURCE_DIR}/common/vec2.cpp
    ${SF_SOURCE_DIR}/common/vec3.cpp
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

# Set library dependencies
set (CPUMSD_LIBS ${MATH_LIB} ${TIME_LIB} ${XML_LIB} ${BOOST_THREAD_LIB} ${NATIVE_THREAD_LIB} ${OPENGL_LIBRARY})

# Set name of the library
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (CPUMSD_LIB "CpuMsd-debug")
else ()
    set (CPUMSD_LIB "CpuMsd")
endif ()

# Final set up to build the CpuMSD library
add_library (${CPUMSD_LIB} SHARED ${CPUMSD_SRCS})
target_link_libraries (${CPUMSD_LIB} ${CPUMSD_LIBS})

if (NOT CMAKE_BUILD_TYPE)
  set_target_properties (${CPUMSD_LIB} PROPERTIES COMPILE_FLAGS "-O2 -Wall")
else ()
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set_target_properties (${CPUMSD_LIB} PROPERTIES COMPILE_FLAGS "-DNDEBUG -O4 -Wall")
  else ()
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set_target_properties (${CPUMSD_LIB} PROPERTIES COMPILE_FLAGS "-g -fno-inline -Wall")
    endif ()
  endif ()
endif ()
