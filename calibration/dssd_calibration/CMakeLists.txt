cmake_minimum_required(VERSION 3.18)
project(calibration)

find_package(cxxopts REQUIRED)

# pad-vetoed thin DSSDs U1, U2, U3, U6 based on protons from 25Si
add_executable(makePadVetoedHists makePadVetoedHists.cpp)
target_include_directories(makePadVetoedHists PRIVATE include)
target_link_libraries(makePadVetoedHists PRIVATE ROOT::Core ROOT::TreePlayer)

add_executable(analysePadVetoedHists analysePadVetoedHists.cpp)
target_include_directories(analysePadVetoedHists PRIVATE include)
target_link_libraries(analysePadVetoedHists PRIVATE ROOT::Core ROOT::TreePlayer ROOT::Spectrum)




# alpha cal. runs, all dssds
add_executable(make3aHists make3aHists.cpp)
target_include_directories(make3aHists PRIVATE include)
target_link_libraries(make3aHists PRIVATE ROOT::Core ROOT::TreePlayer)

add_executable(analyse3aHists analyse3aHists.cpp)
target_include_directories(analyse3aHists PRIVATE include)
target_link_libraries(analyse3aHists PRIVATE ROOT::Core ROOT::TreePlayer ROOT::Spectrum)



# thick DSSDs U4, U5 based on protons from 25Si
add_executable(makeProtonHists makeProtonHists.cpp)
target_include_directories(makeProtonHists PRIVATE include)
target_link_libraries(makeProtonHists PRIVATE ROOT::Core ROOT::TreePlayer)

add_executable(analyseProtonHists analyseProtonHists.cpp)
target_include_directories(analyseProtonHists PRIVATE include)
target_link_libraries(analyseProtonHists PRIVATE ROOT::Core ROOT::TreePlayer ROOT::Spectrum)
