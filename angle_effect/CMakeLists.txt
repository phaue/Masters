cmake_minimum_required(VERSION 3.18)
project(angleeffect)

link_libraries(AUSA::AUSA config++ telescope::telescope GSL::gsl)
file(GLOB HEADERS CONFIGURE_DEPENDS "include/*.h")

add_executable(anglesim angle_effect.cpp ${HEADERS})
add_executable(EnergyEstimation EnergyEstimation.cpp ${HEADERS})
add_executable(PadEE PadEE.cpp ${HEADERS})
add_executable(SolidAngle SolidAngle.cpp ${HEADERS})
