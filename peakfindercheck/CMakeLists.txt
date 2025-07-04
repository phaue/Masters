cmake_minimum_required(VERSION 3.15)
project(peakcheck)

link_libraries(AUSA::AUSA config++ telescope::telescope GSL::gsl)

file(GLOB HEADERS CONFIGURE_DEPENDS "include/*.h")
add_executable(peakcheck peakfinder.cpp ${HEADERS})



target_include_directories(peakcheck PRIVATE
        ${CMAKE_SOURCE_DIR}/cmake           # Include the cmake directory for projectutil.h
)