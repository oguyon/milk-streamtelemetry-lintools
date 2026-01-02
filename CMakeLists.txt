cmake_minimum_required(VERSION 3.10)
project(milk-streamtelemetry-cca C)

find_package(PkgConfig REQUIRED)

pkg_check_modules(CFITSIO REQUIRED cfitsio)
pkg_check_modules(OPENBLAS REQUIRED openblas)
pkg_check_modules(LAPACKE REQUIRED lapacke)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall")

# Common utilities
add_library(common STATIC common.c)
target_include_directories(common PRIVATE ${CFITSIO_INCLUDE_DIRS})
target_link_libraries(common PRIVATE ${CFITSIO_LIBRARIES})

# CCA Program
add_executable(milk-streamtelemetry-cca cca.c)
target_include_directories(milk-streamtelemetry-cca PRIVATE
    ${CFITSIO_INCLUDE_DIRS}
    ${OPENBLAS_INCLUDE_DIRS}
    ${LAPACKE_INCLUDE_DIRS}
)
target_link_libraries(milk-streamtelemetry-cca PRIVATE
    common
    ${CFITSIO_LIBRARIES}
    ${OPENBLAS_LIBRARIES}
    ${LAPACKE_LIBRARIES}
    m
)

# PCA Program
add_executable(milk-streamtelemetry-pca pca.c)
target_include_directories(milk-streamtelemetry-pca PRIVATE
    ${CFITSIO_INCLUDE_DIRS}
    ${OPENBLAS_INCLUDE_DIRS}
    ${LAPACKE_INCLUDE_DIRS}
)
target_link_libraries(milk-streamtelemetry-pca PRIVATE
    common
    ${CFITSIO_LIBRARIES}
    ${OPENBLAS_LIBRARIES}
    ${LAPACKE_LIBRARIES}
    m
)

# Reconstruction Program
add_executable(milk-streamtelemetry-recon recon.c)
target_include_directories(milk-streamtelemetry-recon PRIVATE
    ${CFITSIO_INCLUDE_DIRS}
    ${OPENBLAS_INCLUDE_DIRS}
)
target_link_libraries(milk-streamtelemetry-recon PRIVATE
    common
    ${CFITSIO_LIBRARIES}
    ${OPENBLAS_LIBRARIES}
    m
)
