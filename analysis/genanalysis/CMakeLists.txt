cmake_minimum_required(VERSION 3.18)
project(GeneralAnalysis)

link_libraries(AUSA::AUSA config++ telescope::telescope GSL::gsl)

file(GLOB HEADERS CONFIGURE_DEPENDS "include/*.h")
add_executable(mainanalyzer mainanalyzer.cpp ${HEADERS})
target_include_directories(mainanalyzer PRIVATE
        ${CMAKE_SOURCE_DIR}/cmake           # Include the cmake directory for projectutil.h
        ${CMAKE_SOURCE_DIR}/analysis/genanalysis/include  # Local includes for genanalysis
)