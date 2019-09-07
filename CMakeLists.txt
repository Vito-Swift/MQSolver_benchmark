# Filename:
#       CMakeLists.txt
# Author:
#       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
#
# Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
#
# Description:
#       Primal CMakeLists to build up the project.
#
# Remark:
#       This code is licensed under GNU license.
#       see LICENSE.txt for details.

CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
PROJECT(MQSolver)

SET(CMAKE_CXX_STANDARD_REQUIRED 14)

# Check if avx512 is available on the target machine
INCLUDE(cmake/CHECK_AVX512_PRES.cmake)
IF (AVX512_FOUND)
    MESSAGE(STATUS "[CMakeLists.txt]\tBuilding MQSolver using AVX512 intrinsic set.")
    SET(CMAKE_CXX_FLAGS "-march=native -Ofast -funroll-loops -flto")
    ADD_COMPILE_DEFINITIONS(MODE_AVX512)
    FILE(GLOB MQSOURCE src/mqsolver_avx512/*.cpp)
    FILE(GLOB MQHEADER include/mqsolver_avx512/*.h include/utils.h include/param.h include/mqsolver.h)
ELSE (NOT AVX512_FOUND)
    MESSAGE(STATUS "[CMakeLists.txt]\tBuilding MQSolver using AVX2 intrinsic set.")
    SET(CMAKE_CXX_FLAGS "-march=core-avx2 -Ofast -funroll-loops -flto")
    ADD_COMPILE_DEFINITIONS(MODE_AVX2)
    FILE(GLOB MQSOURCE src/mqsolver_avx2/*.cpp)
    FILE(GLOB MQHEADER include/mqsolver_avx2/*.h include/utils.h include/param.h include/mqsolver.h)
ENDIF (AVX512_FOUND)

# Build library
INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(src)
ADD_LIBRARY(MQLIB ${MQHEADER} ${MQSOURCE})

# Build executables
ADD_EXECUTABLE(mqsolver exec/mqsolver ${MQHEADER} ${MQSOURCE})
TARGET_LINK_LIBRARIES(mqsolver ${MQLIB})

# Copy testing file into build directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/mq46.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/mq46_resident.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
