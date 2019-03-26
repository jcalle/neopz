########## next target ##############

file(GLOB headers *.h)
file(GLOB sources *.cpp)

configure_file(meshes/wellMesh.geo wellMesh.geo COPYONLY)
configure_file(meshes/wellMeshPML.geo wellMeshPML.geo COPYONLY)
configure_file(meshes/concentricMesh.geo concentricMesh.geo COPYONLY)
configure_file(meshes/realWell.geo realWell.geo COPYONLY)
configure_file(meshes/realWellW.geo realWellW.geo COPYONLY)
configure_file(meshes/realWellWFilm.geo realWellWFilm.geo COPYONLY)
configure_file(meshes/axiSymmetricMesh.geo axiSymmetricMesh.geo COPYONLY)
configure_file(meshes/axiRealWell.geo axiRealWell.geo COPYONLY)
configure_file(meshes/axiRealWellW.geo axiRealWellW.geo COPYONLY)
configure_file(meshes/axiRealWellWFilm.geo axiRealWellWFilm.geo COPYONLY)
configure_file(createDirs.sh createDirs.sh COPYONLY)
set(main_SRCS ${sources} ${headers})
add_executable(AcousticsFourierAnalysis ${main_SRCS})
target_link_libraries(AcousticsFourierAnalysis pz)
