# Install script for directory: /gpfs/home/pgersberg/libeigen/Eigen

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/gpfs/home/pgersberg/libeigen/Eigen/QtAlignedMalloc"
    "/gpfs/home/pgersberg/libeigen/Eigen/CholmodSupport"
    "/gpfs/home/pgersberg/libeigen/Eigen/MetisSupport"
    "/gpfs/home/pgersberg/libeigen/Eigen/SuperLUSupport"
    "/gpfs/home/pgersberg/libeigen/Eigen/SPQRSupport"
    "/gpfs/home/pgersberg/libeigen/Eigen/Eigen"
    "/gpfs/home/pgersberg/libeigen/Eigen/Core"
    "/gpfs/home/pgersberg/libeigen/Eigen/StdDeque"
    "/gpfs/home/pgersberg/libeigen/Eigen/SparseQR"
    "/gpfs/home/pgersberg/libeigen/Eigen/SparseCore"
    "/gpfs/home/pgersberg/libeigen/Eigen/Eigenvalues"
    "/gpfs/home/pgersberg/libeigen/Eigen/Geometry"
    "/gpfs/home/pgersberg/libeigen/Eigen/SVD"
    "/gpfs/home/pgersberg/libeigen/Eigen/OrderingMethods"
    "/gpfs/home/pgersberg/libeigen/Eigen/UmfPackSupport"
    "/gpfs/home/pgersberg/libeigen/Eigen/Cholesky"
    "/gpfs/home/pgersberg/libeigen/Eigen/SparseCholesky"
    "/gpfs/home/pgersberg/libeigen/Eigen/StdList"
    "/gpfs/home/pgersberg/libeigen/Eigen/QR"
    "/gpfs/home/pgersberg/libeigen/Eigen/PaStiXSupport"
    "/gpfs/home/pgersberg/libeigen/Eigen/Sparse"
    "/gpfs/home/pgersberg/libeigen/Eigen/Householder"
    "/gpfs/home/pgersberg/libeigen/Eigen/PardisoSupport"
    "/gpfs/home/pgersberg/libeigen/Eigen/SparseLU"
    "/gpfs/home/pgersberg/libeigen/Eigen/StdVector"
    "/gpfs/home/pgersberg/libeigen/Eigen/Jacobi"
    "/gpfs/home/pgersberg/libeigen/Eigen/Dense"
    "/gpfs/home/pgersberg/libeigen/Eigen/IterativeLinearSolvers"
    "/gpfs/home/pgersberg/libeigen/Eigen/LU"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/gpfs/home/pgersberg/libeigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

