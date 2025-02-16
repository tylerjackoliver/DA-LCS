project(ACROBAT LANGUAGES CXX C)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set(Eigen_DIR /usr/local/share/eigen3/cmake)

find_package (Boost REQUIRED)
find_package (dace 2.0.0 REQUIRED)
find_package (Eigen3 REQUIRED)
find_package (MPI REQUIRED)
include_directories( ${Boost_INCLUDE_DIR} )
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_DIR})
add_executable(ACROBAT
	src/main.cpp
)

# Set flags
target_compile_options(ACROBAT PRIVATE -fopenmp -O2 -Wall -std=c++20 -gdwarf-2 -gstrict-dwarf)

# Add linkers
target_include_directories(ACROBAT PRIVATE ${CMAKE_BINARY_DIR}/../include/spice)
target_include_directories(ACROBAT PRIVATE ${Boost_INCLUDE_DIRS})
target_link_directories(ACROBAT PUBLIC ${CMAKE_BINARY_DIR}/../lib)
target_link_libraries(ACROBAT PUBLIC dace::dace_s MPI::MPI_CXX)
target_link_libraries(ACROBAT PRIVATE -fopenmp cspice)
