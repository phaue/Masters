cmake_minimum_required(VERSION 3.15)
project(Masters)

find_package(GSL REQUIRED)
find_package(ROOT REQUIRED)
find_package(AUSALIB REQUIRED)
find_package(telescope REQUIRED)
find_package(libconfig++ REQUIRED)


# Set a global include directory for all subprojects
include_directories(${CMAKE_SOURCE_DIR}/cmake)
include_directories(include)


# Add each analysis subdirectory
add_subdirectory(analysis/genanalysis)
add_subdirectory(analysis/ProtonAnalysis)
add_subdirectory(calibration/pad_calibration)
add_subdirectory(calibration/dssd_calibration)
add_subdirectory(angle_effect)

#This makes it so i can build from the Masters, and simply build all of the analyses at once
# seems quite convenient, however i can still do it the old way inside the specific analysis
