ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (SubStruct SubStruct)
ADD_TEST (SubStruct_Test SubStruct_Test)

#Checking generated files
SET (filename "Cube.vtk")
SET (expectedMD5 "64aa9b689c442c5bf3b0a2b428fc0985")
add_test (SubStruct_Test_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

SET (filename "dohrmann_elastic.scal_vec.0.vtk")
SET (expectedMD5 "faf74fdba1075ee2bc3b4917bfe7eb5e")
add_test (SubStruct_Test_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

SET (filename "malhaPZ1BC.txt")
SET (expectedMD5 "fdc16e681800611f9730a5ffcf9d231b")
add_test (SubStruct_Test_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#Checking generated files
#SET (filename "CheckPoint1.txt")
#SET (expectedMD5 "edd1c50d9a67f266ae8f75f04e36d474")
#add_test (SubStruct_CheckPoint1.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "CheckPoint2.txt")
#SET (expectedMD5 "33469e4023d4cee6042a872654e0487d")
#add_test (SubStruct_CheckPoint2.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "CheckPoint3.txt")
#SET (expectedMD5 "a94515682eaec56b21414b4cc5f994d1")
#add_test (SubStruct_CheckPoint3.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "CoarseMatrix.vtk")
#SET (expectedMD5 "9ebe6c8f7c59803b4e090823372560ad")
#add_test (SubStruct_CoarseMatrix.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "Cube.vtk")
#SET (expectedMD5 "64aa9b689c442c5bf3b0a2b428fc0985")
#add_test (SubStruct_Cube.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "dohrmann_visco.scal_vec.0.vtk")
#SET (expectedMD5 "5e94a14e6871df404ce792937b548d72")
#add_test (SubStruct_dohrmann_visco.scal_vec.0.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "malhaPZ1BC.txt")
#SET (expectedMD5 "fdc16e681800611f9730a5ffcf9d231b")
#add_test (SubStruct_malhaPZ1BC.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "partition.vtk")
#SET (expectedMD5 "f0212d1b5532f48d8ebbb45a4366a67c")
#add_test (SubStruct_partition.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "partitionbefore.vtk")
#SET (expectedMD5 "5c856bcaab0bde24e5f4ace219c36857")
#add_test (SubStruct_partitionbefore.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "PointMesh.vtk")
#SET (expectedMD5 "d9f2bf12bfaf18e0f49a32736d6963ed")
#add_test (SubStruct_PointMesh.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrix79.vtk")
#SET (expectedMD5 "442d28ac8e338ee416d2cb5bad3c9f13")
#add_test (SubStruct_SubMatrix79.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrix80.vtk")
#SET (expectedMD5 "ab6fa7c83836c1851dd962d4b5de7a08")
#add_test (SubStruct_SubMatrix80.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrix81.vtk")
#SET (expectedMD5 "37815aabd7333c1b0ee8f8f7aedbd161")
#add_test (SubStruct_SubMatrix81.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrix82.vtk")
#SET (expectedMD5 "8f1d28ef2f1e8c9ad60d187b6d5a3c67")
#add_test (SubStruct_SubMatrix82.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrix83.vtk")
#SET (expectedMD5 "5dc8e911767abbb0f23765100116e454")
#add_test (SubStruct_SubMatrix83.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrix84.vtk")
#SET (expectedMD5 "e5a9aabc517c05f6aeac9c9d9a3a04c8")
#add_test (SubStruct_SubMatrix84.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrixInternal79.vtk")
#SET (expectedMD5 "5da92605d0bf66aa048fdd628366fbf3")
#add_test (SubStruct_SubMatrixInternal79.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrixInternal80.vtk")
#SET (expectedMD5 "ab6fa7c83836c1851dd962d4b5de7a08")
#add_test (SubStruct_SubMatrixInternal80.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrixInternal81.vtk")
#SET (expectedMD5 "75c393a07a6ed6b6c22e81c8afd27b20")
#add_test (SubStruct_SubMatrixInternal81.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrixInternal82.vtk")
#SET (expectedMD5 "8f1d28ef2f1e8c9ad60d187b6d5a3c67")
#add_test (SubStruct_SubMatrixInternal82.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrixInternal83.vtk")
#SET (expectedMD5 "d0d7f1e2c1c6edbae07af0b37a6b722e")
#add_test (SubStruct_SubMatrixInternal83.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "SubMatrixInternal84.vtk")
#SET (expectedMD5 "7b06fff5a06321d3abb2aee6b94ab0ed")
#add_test (SubStruct_SubMatrixInternal84.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
