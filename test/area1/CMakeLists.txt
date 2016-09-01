# CMake file for test-area1
project (SF_TEST1 CXX)

# Set source file names
set (TEST1_SRCS
    ${SF_SOURCE_DIR}/common/vec3.cpp
    ${SF_SOURCE_DIR}/common/vec4.cpp
    src/EM_common.cpp
    src/EM_FEMSubmesh.cpp
    src/EM_FEMMesh.cpp
    src/EM_MSDMesh.cpp
    src/main.cpp)

include_directories (${SF_SOURCE_DIR}/common ./inc)

set (TEST1_LIBS ${MATH_LIB})

# Set final name of executable
add_executable (test1 ${TEST1_SRCS})
target_link_libraries (test1 ${TEST1_LIBS})