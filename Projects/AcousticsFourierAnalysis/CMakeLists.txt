########## next target ##############

file(GLOB headers *.h)
file(GLOB sources *.cpp)

configure_file(wellMesh.geo wellMesh.geo COPYONLY)
configure_file(wellMeshPML.geo wellMeshPML.geo COPYONLY)
configure_file(concentricMesh.geo concentricMesh.geo COPYONLY)
configure_file(realWell.geo realWell.geo COPYONLY)
set(main_SRCS ${sources} ${headers})
add_executable(AcousticsFourierAnalysis ${main_SRCS})
target_link_libraries(AcousticsFourierAnalysis pz)
