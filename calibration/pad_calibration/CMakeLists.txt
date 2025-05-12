cmake_minimum_required(VERSION 3.18)
project(PadCalibration)

link_libraries(AUSA::AUSA config++ telescope::telescope GSL::gsl)

file(GLOB HEADERS CONFIGURE_DEPENDS "include/*.h")
add_executable(padcal PadCalibration.cpp ${HEADERS})
add_executable(fitter fitter.cpp ${HEADERS})
target_include_directories(padcal PRIVATE
        ${CMAKE_SOURCE_DIR}/cmake           # Include the cmake directory for projectutil.h
        ${CMAKE_SOURCE_DIR}/calibration/pad_calibration/include  # Local includes for genanalysis
)
target_include_directories(fitter PRIVATE
        ${CMAKE_SOURCE_DIR}/cmake           # Include the cmake directory for projectutil.h
        ${CMAKE_SOURCE_DIR}/calibration/pad_calibration/include  # Local includes for genanalysis
)