# CMake file for test-area1
project (SF_TEST2 CXX)

# Set source file names
set (TEST2_SRCS
    ${SF_SOURCE_DIR}/common/Collide/lineTriCollide.cpp
    ${SF_SOURCE_DIR}/common/Collide/triTriCollide.cpp
    src/main.cpp)

include_directories (${SF_SOURCE_DIR}/common ./inc)

set (TEST2_LIBS ${BOOST_THREAD_LIB} ${MATH_LIB})

# Set final name of executable
add_executable (test2 ${TEST2_SRCS})
target_link_libraries (test2 ${TEST2_LIBS})
