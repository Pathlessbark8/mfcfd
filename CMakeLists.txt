cmake_minimum_required(VERSION 3.10)
project(mfcfd LANGUAGES Fortran)
ENABLE_LANGUAGE(Fortran)
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_compile_options(-march=native -O3 -ffree-line-length-512)

#------------------------------------------------------------------------------#
# PETSc support
#------------------------------------------------------------------------------#
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
set(LIBS ${LIBS} ${PETSC_LIBRARIES})

#------------------------------------------------------------------------------#
# HDF5 support
#------------------------------------------------------------------------------#
find_package(HDF5 COMPONENTS C CXX Fortran REQUIRED)
link_directories( ${HDF5_LIBRARY_DIRS})
include_directories(${HDF5_Fortran_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_Fortran_LIBRARIES})

if( MFCFD EQUAL SERIAL )
message("Building Serial MeshFree Solver")
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_serial/*.F90")
add_executable(execname ${MFCFD_SOURCES})
target_link_libraries(execname ${LIBS})

else()
message("Building Serialsd MeshFree Solver")
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_serial/*.F90")
add_executable(execname ${MFCFD_SOURCES})
target_link_libraries(execname ${LIBS})

endif()