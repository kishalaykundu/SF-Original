# CMake file for Engine/Driver
project (DRIVER CXX)

if (NOT WITH_VECTOR3 OR WITH_VECTOR3 STREQUAL "OFF")
    set (DRIVER_SRCS ${SF_SOURCE_DIR}/common/vec4.cpp)
else ()
    set (DRIVER_SRCS ${SF_SOURCE_DIR}/common/vec3.cpp)
endif ()

set (DRIVER_SRCS ${DRIVER_SRCS}
    ${SF_GRAPHICS_DIR}/Display/src/Display.cpp
    src/Resource.cpp
    src/Driver.cpp
    src/main.cpp)

# Set include directory paths
include_directories (./inc ${SF_SOURCE_DIR}/common ${SF_GRAPHICS_DIR}/Display/inc/)

# platform specific include paths
if (APPLE)
	find_library (PNG_LIB png12 /opt/local/lib/ REQUIRED)
    include_directories (/opt/local/include/)
else ()
    if (WIN32)

    else ()
    	find_library (PNG_LIB png12 REQUIRED)
        include_directories (/usr/include/ /usr/include/libxml2)
    endif ()
endif ()

# find Display system
set (DRIVER_REQUIRED_LIBS ${DRIVER_REQUIRED_LIBS} ${OPENGL_LIBRARY} ${GLU_LIBRARY} ${GLUT_LIB} ${PNG_LIB})

# set essential library links
set (DRIVER_REQUIRED_LIBS ${DRIVER_REQUIRED_LIBS} ${XML_LIB} ${MATH_LIB} ${BOOST_THREAD_LIB} ${NATIVE_THREAD_LIB})

# set and link target
add_executable (runsim ${DRIVER_SRCS})
target_link_libraries (runsim ${DRIVER_REQUIRED_LIBS})

if (NOT CMAKE_BUILD_TYPE)
  set_target_properties (runsim PROPERTIES COMPILE_FLAGS "-O2 -Wall")
else ()
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set_target_properties (runsim PROPERTIES COMPILE_FLAGS "-DNDEBUG -O4 -Wall")
  else ()
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set_target_properties (runsim PROPERTIES COMPILE_FLAGS "-g -fno-inline -Wall")
    endif ()
  endif ()
endif ()
