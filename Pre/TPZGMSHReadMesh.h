/**
 * classroom that manages the manipulation of geometric meshes 
 * generated by the GMSH in order to be used for the NeoPZ
 * the archive generated for the GMSH must contain volume elements 
 * and contour elements,
 * the data of these elements must be:
 * - number of nodes 
 * - co-ordinated of the nodes 
 * - number of elements 
 * - incidences of the elements
 */

#ifndef TPZGMSH_HPP
#define TPZGMSH_HPP

class TPZGeoMesh;
class TPZGeoEl;
class TPZGeoElSide;
#include "pzstack.h"
#include <fstream>
/*************/

class TPZGMSHReadMesh {

  /**
   * archive generated for the GMSH to be used (interpreted) inside of the NeoPZ
   */
	std::ifstream fInGSMHGeoMesh;

  /**
   * defined geometric mesh in the NeoPZ
   */
  TPZGeoMesh *fGeoMesh;

 public:
  
  TPZGMSHReadMesh(TPZGeoMesh *gmesh);
  //TPZGMSHReadMesh(char *meshGMSH,TPZGeoMesh *gmesh);
  
  ~TPZGMSHReadMesh(){};

  /*
   * readings of meshes 2D, the contour and volume elements are returned in  
   * elemlist and elembclist respectively 
   * for a rectangular mesh the CC are:
   * 1: for it wing
   * 2: left lateral contour 
   * 3: inferior contour 
   * 4: right lateral contour
   * 5: superior contour
   * the plan is with axle X for right and axle Y for top
   */
  void ReadMesh2D(char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist);
  void ReadMesh2D2(char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist);
  /*
   * readings of meshes 2D, the contour and volume elements are returned in  
   * elemlist and elembclist respectively 
   * for one mesh of tetrahedrons the CC are
   * 1: for it wing
   * 2: left lateral contour 
   * 3: inferior contour 
   * 4: right lateral contour
   * 5: superior contour
   * 6: posterior contour (Z = 0) 
   * 7: previous contour (Z > 0)
   * taking as base that the domain is with axle X for right, axle Y for top and axle -Z for the deep one
   */
  void ReadMesh3D(char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist);

  /**
   * it rearranges the nodal numeration given by the GMSH of sequential form
   */
  void Resequence(TPZStack<int> &Indexes,char *meshfile);

  /*
   * it prints in the exit defined for out the characteristics of the geometric mesh 
   * created by the NeoPZ based on the mesh generated for the GMSH
   */
	void PrintGeoMesh(std::ostream &out = std::cout);
  
};

#endif
