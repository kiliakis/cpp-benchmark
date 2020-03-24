# Requirements list
cmake_minimum_required(VERSION 3.0.2)
project(RandomGeneratorsBenchmark CXX)
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

#set(CMAKE_BUILD_TYPE Release)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
      set(CMAKE_CXX_FLAGS         "-std=c++11 -fopenmp -Wall -pedantic -fno-rtti -Wunused -fopt-info-vec=report.txt")
      set(CMAKE_CXX_FLAGS_DEBUG   "-Ofast -g3") #-fsanitize=address")
      set(CMAKE_CXX_FLAGS_RELEASE "-march=native")
      
      CHECK_CXX_COMPILER_FLAG(-Ofast COMPILER_HAS_OFAST)
      if(COMPILER_HAS_OFAST)
         set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
      else()
         set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      endif()
      
      set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")
endif()

include_directories(include /usr/local/include)
# include_directories($ENV{HOME}/opt/include)  
# link_directories(/usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
# link_directories("$ENV{HOME}/opt/lib")

find_package(OpenMP REQUIRED)
set(LIBRARIES m)


file(GLOB_RECURSE SOURCES
        "src/*"
        "include/*"
        )

set(SOURCE_FILES ${SOURCES})

# foreach(source ${SOURCE_FILES})
#   set_source_files_properties(${source} PROPERTIES COMPILE_FLAGS "-fopt-info-vec=report.txt")
# endforeach()

add_library(blond SHARED ${SOURCE_FILES})
target_compile_definitions(blond PRIVATE SHARED=1)
target_link_libraries(blond ${LIBRARIES})
LIST(APPEND LIBRARIES blond)




file(GLOB EXECUTABLE_SOURCES "tests/*.cpp")
set(EXECUTABLE_FILES ${EXECUTABLE_SOURCES})

foreach(src ${EXECUTABLE_FILES})
   file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
   string( REPLACE ".cpp" "" testname ${testfile})
   string( REPLACE "tests/" "" testname ${testname})
   add_executable( ${testname} ${src})
   target_link_libraries( ${testname} ${LIBRARIES})
endforeach()

