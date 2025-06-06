cmake_minimum_required(VERSION 3.18)
project(ProtonAnalysis)

link_libraries(AUSA::AUSA config++ telescope::telescope GSL::gsl)

file(GLOB HEADERS CONFIGURE_DEPENDS "include/*.h")
add_executable(makeHists makeHists.cpp ${HEADERS})
target_include_directories(makeHists PRIVATE
        ${CMAKE_SOURCE_DIR}/cmake           # Include the cmake directory for projectutil.h
        ${CMAKE_SOURCE_DIR}/analysis/ProtonAnalysis/include  # Local includes for genanalysis
)