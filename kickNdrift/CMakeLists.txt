# Requirements list
cmake_minimum_required(VERSION 3.0.2)
project(kickNDriftBench CXX)
include(CheckCXXCompilerFlag)

# TODO is 4.8 really the minimum requirement?
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0")
    message(FATAL_ERROR "ERROR: Insufficient compiler version, need at lest 4.8 to continue")
endif()

check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CXX11)
if(!COMPILER_HAS_CXX11)
   message(FATAL_ERROR "ERROR: Compiler is lacking support for c++11")
endif()



if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_FLAGS "-std=gnu++11 -fopenmp -DUSE_OMP -Wall -pedantic -fno-rtti -Wunused -g")
set(CMAKE_CXX_FLAGS_DEBUG   "-g3") #-fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-march=native")
# set(LIBRARIES m)

CHECK_CXX_COMPILER_FLAG(-Ofast COMPILER_HAS_OFAST)
if(COMPILER_HAS_OFAST)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -ffast-math")
else()
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

#set(CMAKE_BUILD_TYPE Release)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
   set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
   # set(LIBRARIES "mkl_intel_ilp64" "mkl_intel_thread" "mkl_core" "iomp5" "pthread" "m" "dl")
   # set(LIBRARIES "mkl_intel_ilp64" "mkl_tbb_thread" "mkl_core" "tbb" "stdc++" "pthread" "m" "dl")
endif()

# Magic to set GCC-specific compile flags (to turn on optimisation).
#set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -msse3 -mavx")
#add_definitions( -DSSE3 )
#add_definitions( -DAVX )

# set(LIBRARIES m mkl_)

# include_directories(include /afs/cern.ch/sw/IntelSoftware/linux/x86_64/xe2017/compilers_and_libraries/linux/tbb/include/tbb)
include_directories(include)


file(GLOB EXECUTABLE_SOURCES "src/*.cpp")
set(EXECUTABLE_FILES ${EXECUTABLE_SOURCES})

foreach(src ${EXECUTABLE_FILES})
   file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
   string( REPLACE ".cpp" "" testname ${testfile})
   string( REPLACE "src/" "" testname ${testname})
   add_executable( ${testname} ${src})
   target_link_libraries( ${testname} ${LIBRARIES})
endforeach()

