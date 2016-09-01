# Top level file for the CU_MSD (CUDA accelerated eXtended Finite Element Method) module
project (SF_RIGID CXX)

if (NOT WITH_VECTOR3 OR WITH_VECTOR3 STREQUAL "OFF")
    set (RIGID_SRCS ${SF_SOURCE_DIR}/common/vec4.cpp)
else ()
    set (RIGID_SRCS ${SF_SOURCE_DIR}/common/vec3.cpp)
endif ()

set (RIGID_SRCS ${RIGID_SRCS}
    ${SF_SOURCE_DIR}/common/vec2.cpp
    ${SF_SOURCE_DIR}/common/mat4x4.cpp
    ${SF_SOURCE_DIR}/common/GL/common.cpp
    src/Common.cpp
    src/Mesh.cpp
    src/Plugin.cpp)

# Set platform-independent include paths
include_directories ( ./inc
    ${SF_SOURCE_DIR}/common/
    ${SF_ENGINE_DIR}/inc/
    ${SF_GRAPHICS_DIR}/Display/inc/
    /usr/include/libxml2/)

# Set library dependencies
set (RIGID_LIBS ${MATH_LIB} ${XML_LIB} ${BOOST_THREAD_LIB} ${NATIVE_THREAD_LIB} ${OPENGL_LIBRARY})

# Set name of the library
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (RIGID_LIB "Rigid-debug")
else ()
    set (RIGID_LIB "Rigid")
endif ()

# Final set up to build the Rigid library
add_library (${RIGID_LIB} SHARED ${RIGID_SRCS})
target_link_libraries (${RIGID_LIB} ${RIGID_LIBS})

if (NOT CMAKE_BUILD_TYPE)
  set_target_properties (${RIGID_LIB} PROPERTIES COMPILE_FLAGS "-O2 -Wall")
else ()
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set_target_properties (${RIGID_LIB} PROPERTIES COMPILE_FLAGS "-DNDEBUG -O4 -Wall")
  else ()
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set_target_properties (${RIGID_LIB} PROPERTIES COMPILE_FLAGS "-g -fno-inline -Wall")
    endif ()
  endif ()
endif ()
