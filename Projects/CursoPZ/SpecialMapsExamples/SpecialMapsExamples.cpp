#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
//special maps
//#include "TPZCylinder.h"//@TODO: check if this class is still useful
#include "TPZQuadTorus.h"
#include "TPZTriangleTorus.h"
#include "tpzellipse3d.h"
#include "tpzquadraticline.h"
#include "tpzquadraticpyramid.h"
#include "tpzquadratictetra.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"
#include "TPZWavyLine.h"
#include "tpzarc3d.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"
//#include "tpzblendnaca.h"//@TODO: create InsertExampleElement for this mapping

template <class TGeo>
void CreateSampleElement(const int &nref);

int main(int argc, char *argv[]) {
    const int nRef = 4;
    CreateSampleElement<pzgeom::TPZQuadTorus>(nRef);
//    CreateSampleElement<pzgeom::TPZArc3D>(nRef);
}

template <class TGeo>
void CreateSampleElement(const int &nref){
    std::string elName = TGeo::TypeName();

    std::cout<<"Starting mesh refinement of element type: "<<elName<<std::endl;

    const int matId = 1;
    TPZManVector<REAL,3> lowerCorner(3,0);
    TPZManVector<REAL,3> size(3,1);
    TPZGeoMesh gmesh;
    TGeo::InsertExampleElement(gmesh, matId, lowerCorner, size);
    gmesh.BuildConnectivity();

    {
        TPZVec<TPZGeoEl *> sons;
        for(int iDiv = 0; iDiv < nref; iDiv++){
            const int nel = gmesh.NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *geo = gmesh.Element(iel);
                if (geo->NSubElements() == 0) {
                    geo->Divide(sons);
                }
            }
        }
    }
    std::ofstream fileName(elName+".vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&gmesh, fileName);
}