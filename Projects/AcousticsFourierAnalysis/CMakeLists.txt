########## next target ##############

file(GLOB headers *.h)
file(GLOB sources *.cpp)

configure_file(wellMesh.geo wellMesh.geo COPYONLY)
set(main_SRCS ${sources} ${headers})
add_executable(AcousticsFourierAnalysis ${main_SRCS})
target_link_libraries(AcousticsFourierAnalysis pz)
