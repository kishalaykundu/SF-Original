# Top Level CMakeLists.txt file for the Integrate system
project (SIMULATE)

# Minimum cmake version set to 2.4.5 (using version 2.6.4)
cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

# Allows non-verbose formats for cmake loop constructs
set (CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if (COMMAND cmake_policy)
    cmake_policy (SET CMP0003 NEW)
endif (COMMAND cmake_policy)

# Set shared libraries to be built
set (BUILD_SHARED_LIBS ON)

######################### SET UP SF DIRECTORY STRUCTURES #########################

# Set some cmake variables to be used
set (SF_SOURCE_DIR ${CMAKE_SOURCE_DIR})
set (SF_ENGINE_DIR ${CMAKE_SOURCE_DIR}/driver/)
set (SF_PLUGINS_DIR ${CMAKE_SOURCE_DIR}/plugins/)
set (SF_TEST_DIR ${CMAKE_SOURCE_DIR}/test/)
set (SF_UTILITIES_DIR ${CMAKE_SOURCE_DIR}/utils/)

# Set default directories for haptics (to be precompiled with the engine)
set (SF_GRAPHICS_DIR ${CMAKE_SOURCE_DIR}/plugins/graphics/ogl)
set (SF_HAPTICS_DIR ${CMAKE_SOURCE_DIR}/plugins/haptics/)

#######################################################################################


###################### IDENTIFY LIBRARIES REQUIRED BY SF ######################

# Common paths and libraries used through the project
find_library (MATH_LIB m REQUIRED)
find_library (DL_LIB dl REQUIRED)
find_library (TIME_LIB boost_date_time REQUIRED)
find_library (XML_LIB xml2 REQUIRED)
find_library (BOOST_THREAD_LIB boost_thread REQUIRED)

# platform-specific required libraries
if (APPLE) # mac os x
    find_library (NATIVE_THREAD_LIB pthread REQUIRED)
else ()

    if (WIN32) # windows

    else () # linux
    	find_library (NATIVE_THREAD_LIB pthread REQUIRED)
    endif ()

endif ()

# find OpenGL
if (APPLE)
    find_package (OPENGL REQUIRED)
else ()
    find_package (OpenGL REQUIRED)
endif ()

# find GLU
find_package (GLU REQUIRED)

# find GLUT
find_library (GLUT_LIB glut REQUIRED)

mark_as_advanced (FORCE OPENGL_LIBRARY GLU_LIBRARY GLUT_LIB)
mark_as_advanced (MATH_LIB DL_LIB TIME_LIB XML_LIB BOOST_THREAD_LIB NATIVE_THREAD_LIB)

#######################################################################################

# List of all subdirectories to be parsed
add_subdirectory (${SF_ENGINE_DIR})
add_subdirectory (${SF_PLUGINS_DIR})
add_subdirectory (${SF_TEST_DIR})
add_subdirectory (${SF_UTILITIES_DIR})

###################### SET UP CUSTOM FLAGS TO BE PASSED BY CMAKE ######################

# set option flag for double-precision
# to set, call with command -DDOUBLE_PRECISION=ON
if (NOT DOUBLE_PRECISION OR DOUBLE_PRECISION STREQUAL "OFF")
    option (SF_DOUBLE_PRECISION "Disable double precision" OFF)
else ()
    option (SF_DOUBLE_PRECISION "Enable double precision" ON)
endif ()

# set option flag for 3-D vector usage (default: 4-D homogeneous co-ordinates)
# to set, call with command -DWITH_VECTOR3=ON
if (NOT WITH_VECTOR3)
    option (SF_VECTOR3_ENABLED "3D Vector" OFF)
    option (SF_VECTOR4_ENABLED "4D Vector" ON)
else ()
    option (SF_VECTOR3_ENABLED "3D Vector" ON)
    option (SF_VECTOR4_ENABLED "4D Vector" OFF)
endif ()

# set option for print to terminal
# to set, call with command -DNO_PRINT=ON (default: OFF)
if (NOT NO_PRINT OR NO_PRINT STREQUAL "OFF")
    option (SF_NO_PRINT "No print to terminal" OFF)
else ()
    option (SF_NO_PRINT "No print to terminal" ON)
endif ()

# IMPORTANT: This line should be the last line after all option flags are set
configure_file (${SF_SOURCE_DIR}/common/config.h.in ${SF_SOURCE_DIR}/common/config.h)

#######################################################################################

# Installation options
set (CMAKE_INSTALL_PREFIX "/home/kish1/Projects")

#######################################################################################

###################### ECHO CUSTOM CMAKE OPTIONS ######################
message (STATUS "")
message (STATUS "CMAKE OPTIONS FOR SF")
message (STATUS "\t-DCMAKE_BUILD_TYPE=Debug/None/Release (default: None)")
message (STATUS "\t-DDOUBLE_PRECISION=ON/OFF (default: OFF)")
message (STATUS "\t-DWITH_VECTOR3=ON/OFF (default: Vector4)")
message (STATUS "\t-DNO_PRINT=ON/OFF (default: OFF)")
message (STATUS "")
