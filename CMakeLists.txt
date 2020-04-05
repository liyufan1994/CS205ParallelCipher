cmake_minimum_required(VERSION 3.15)
project(Denigma)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)
include_directories(lib)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)

#set (CMAKE_CXX_FLAGS "-fopenmp")
#set (CMAKE_C_FLAGS "-fopenmp")


add_library(ArrayUtils SHARED src/ArrayUtilities.cpp)
add_library(Sampling SHARED src/Randomize.cpp)
link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable(Denigma
        src/CypherUtilities.cpp
        src/decipher.cpp
        src/decipher.h
        src/cipherMCMC.cpp
        src/FineSearch.cpp
        src/ArrayUtilities.h
        src/Randomize.h)


target_link_libraries(Denigma ArrayUtils Sampling)

add_executable(Ising
        src/ArrayUtilities.h
        src/Randomize.h
        src/Ising.cpp
        src/Ising.h
        src/IsingMCMC.cpp)

target_link_libraries(Ising ArrayUtils Sampling)
