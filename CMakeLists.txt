#MFCFD CMAKE SCRIPT

cmake_minimum_required(VERSION 3.10)
project(mfcfd)
ENABLE_LANGUAGE(Fortran)
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if( NOT "${MFCFD}" MATCHES "CUDA" )
#------------------------------------------------------------------------------#
# PETSc support
#------------------------------------------------------------------------------#
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
set(LIBS ${LIBS} ${PETSC_LIBRARIES})

else()
#------------------------------------------------------------------------------#
# CUDA support
#------------------------------------------------------------------------------#
find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CUDA_LIBRARIES})

endif()

#------------------------------------------------------------------------------#
# HDF5 support
#------------------------------------------------------------------------------#
find_package(HDF5 COMPONENTS C CXX Fortran REQUIRED)
link_directories( ${HDF5_LIBRARY_DIRS})
include_directories(${HDF5_Fortran_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_Fortran_LIBRARIES})

if( "${MFCFD}" MATCHES "SERIAL" )
message("Building Serial MeshFree Solver")
add_compile_options(-march=native -O3 -ffree-line-length-512)
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_serial/*.F90")
add_executable(execname ${MFCFD_SOURCES})
target_link_libraries(execname ${LIBS})

elseif( "${MFCFD}" MATCHES "TANGENT" )
message("Building Tangent Meshfree Solver")
add_compile_options(-march=native -O3 -ffree-line-length-512)
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_forward/*.F90")
FILE(GLOB MFCFD_SENSE_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_forward/sensitivity_files/*.F90")
add_executable(execname ${MFCFD_SOURCES} ${MFCFD_SENSE_SOURCES})
target_link_libraries(execname ${LIBS})

elseif( "${MFCFD}" MATCHES "ADJOINT" )
message("Building Adjoint Meshfree Solver")
add_compile_options(-march=native -O3)
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_adjoint/*.F90")
FILE(GLOB MFCFD_SENSE_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_adjoint/sensitivity_files/*.F90")
set(MFCFD_AD_SOURCES ${CMAKE_SOURCE_DIR}/src_mpi_adjoint/ADFirstAidKit/adBuffer.f ${CMAKE_SOURCE_DIR}/src_mpi_adjoint/ADFirstAidKit/adStack.c ${CMAKE_SOURCE_DIR}/src_mpi_adjoint/revolve.f)
#message(${MFCFD_AD_SOURCES})
set_source_files_properties(${MFCFD_SOURCES} PROPERTIES COMPILE_FLAGS -ffree-line-length-512)
set_source_files_properties(${MFCFD_SENSE_SOURCES} PROPERTIES COMPILE_FLAGS -ffree-line-length-512)
set_source_files_properties(${CMAKE_SOURCE_DIR}/src_mpi_adjoint/revolve.f PROPERTIES COMPILE_FLAGS -fno-range-check)
set_source_files_properties(${CMAKE_SOURCE_DIR}/src_mpi_adjoint/ADFirstAidKit/adBuffer.f PROPERTIES COMPILE_FLAGS -fno-range-check)
add_executable(execname ${MFCFD_SOURCES} ${MFCFD_SENSE_SOURCES} ${MFCFD_AD_SOURCES})
target_link_libraries(execname ${LIBS})

elseif( "${MFCFD}" MATCHES "CUDA" )
message("Building Cuda Meshfree Solver")
#add_compile_options(-O3 -Mcuda=nordc -Mcuda=ptxinfo)
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_cuda/*.F90")
FILE(GLOB MFCFD_CUDA_SOURCES "${CMAKE_SOURCE_DIR}/src_cuda/*.cuf")
set(CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Mcuda=nordc")
add_executable(execname ${MFCFD_SOURCES} ${MFCFD_CUDA_SOURCES})
target_link_libraries(execname ${CUDA_cusparse_LIBRARY} ${LIBS})

else()
message("Building Default Serial MeshFree Solver")
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_serial/*.F90")
add_executable(execname ${MFCFD_SOURCES})
target_link_libraries(execname ${LIBS})

endif()
