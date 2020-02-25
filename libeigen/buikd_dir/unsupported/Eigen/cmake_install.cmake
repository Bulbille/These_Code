# Install script for directory: /gpfs/home/pgersberg/libeigen/unsupported/Eigen

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/AdolcForward"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/AlignedVector3"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/ArpackSupport"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/AutoDiff"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/BVH"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/EulerAngles"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/FFT"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/IterativeSolvers"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/KroneckerProduct"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/LevenbergMarquardt"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/MatrixFunctions"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/MoreVectorization"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/MPRealSupport"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/NonLinearOptimization"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/NumericalDiff"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/OpenGLSupport"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/Polynomials"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/Skyline"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/SparseExtra"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/SpecialFunctions"
    "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/Splines"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/gpfs/home/pgersberg/libeigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/gpfs/home/pgersberg/libeigen/buikd_dir/unsupported/Eigen/CXX11/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

