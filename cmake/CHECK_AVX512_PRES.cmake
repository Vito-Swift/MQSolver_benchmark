# Filename:
#       Check_AVX512_PRES.cmake
# Author:
#       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
#
# Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
#
# Description:
#       Check if AVX512 instructions are available on the
#       machine where the project is compiled.
#
# Remark:
#       This code is licensed under GNU license
#       see LICENSE.txt for details

MESSAGE(STATUS "[check_avx512_pres.cmake]\tdetecting system architecture...")

IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    EXEC_PROGRAM(cat ARGS "/proc/cpuinfo" OUTPUT_VARIABLE CPUINFO)

    # Check AVX512 Presence
    STRING(REGEX REPLACE "^.*(avx512).*$" "\\1" SSE_FLAGS ${CPUINFO})
    STRING(COMPARE EQUAL "avx512" "${SSE_FLAGS}" AVX512_TRUE)
    IF (AVX512_TRUE)
        set(AVX512_FOUND true CACHE BOOL "AVX512 available on host")
    ELSE (AVX512_TRUE)
        set(AVX512_FOUND false CACHE BOOL "AVX512 available on host")
    ENDIF (AVX512_TRUE)

    # Check AVX2 Presence
    STRING(REGEX REPLACE "^.*(avx2).*$" "\\1" SSE_FLAGS ${CPUINFO})
    STRING(COMPARE EQUAL "avx2" "${SSE_FLAGS}" AVX2_TRUE)
    IF (AVX2_TRUE)
        set(AVX2_FOUND true CACHE BOOL "AVX2 available on host")
    ELSE (AVX2_TRUE)
        set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
    ENDIF (AVX2_TRUE)
    #TODO: Add more operating system like Windows, MacOS, etc.
ELSE (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(AVX2_FOUND true CACHE BOOL "AVX2 available on host")
    set(AVX512_FOUND false CACHE BOOL "AVX512 available on host")
ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")

IF (NOT AVX512_FOUND)
    MESSAGE(STATUS "[check_avx512_pres.cmake]\tcannot find hardware support for AVX512 on this machine.")
ELSEIF (AVX512_FOUND)
    MESSAGE(STATUS "[check_avx512_pres.cmake]\tfound AVX2 available on the target machine.")
ENDIF (NOT AVX512_FOUND)
IF (NOT AVX2_FOUND)
    MESSAGE(FATAL_ERROR "[check_avx512_pres.cmake]\tcannot find hardware support for AVX2 on this machine.\nAborting...")
ELSEIF (AVX2_FOUND)
    MESSAGE(STATUS "[check_avx512_pres.cmake]\tfound AVX2 available on the target machine.")
ENDIF (NOT AVX2_FOUND)