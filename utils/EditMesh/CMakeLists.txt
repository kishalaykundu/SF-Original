# CMake file for test-area1
project (SF_EDITMESH CXX)

# Set source file names
set (EDITMESH_SRCS
    ${SF_SOURCE_DIR}/common/vec3.cpp
    ${SF_SOURCE_DIR}/common/vec4.cpp
    src/EM_common.cpp
    src/EM_FEMSubmesh.cpp
    src/EM_FEMMesh.cpp
    src/EM_MSDMesh.cpp
    src/main.cpp)

include_directories (${SF_SOURCE_DIR}/common ./inc)

set (EDITMESH_LIBS ${MATH_LIB})

# Set final name of executable
add_executable (edit-mesh ${EDITMESH_SRCS})
target_link_libraries (edit-mesh ${EDITMESH_LIBS})

if (NOT CMAKE_BUILD_TYPE)
  set_target_properties (edit-mesh PROPERTIES COMPILE_FLAGS "-O2 -Wall")
else ()
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set_target_properties (edit-mesh PROPERTIES COMPILE_FLAGS "-DNDEBUG -O4 -Wall")
  endif ()
endif ()
