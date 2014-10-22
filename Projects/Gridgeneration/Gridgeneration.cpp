#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzmatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "time.h"
#include "pzconvectionproblem.h"
#include "pzmultiphase.h"
#include "pzl2projection.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"


#include <time.h>
#include <stdio.h>

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"
#include <cmath>



void Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X);
void Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X);
void Parametricfunction3(const TPZVec<STATE> &par, TPZVec<STATE> &X);

bool ftriang = false;
REAL angle = 0.0*M_PI/4.0;

int main()
{   

    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
    GeoMesh1->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh1->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
    GeoMesh1->BuildConnectivity();
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMesh1.txt");
        GeoMesh1->Print(argument);
        std::ofstream Dummyfile("GeometricMesh1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh1,Dummyfile, true);
    }
    
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh1);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc = new TPZDummyFunction<STATE>(Parametricfunction);
    CreateGridFrom.SetParametricFunction(ParFunc);
    REAL t=0.0;
    REAL dt = 0.05;
    int n = 50;
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh2 = CreateGridFrom.ComputeExtrusion(t, dt, n);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew2.txt");
        GeoMesh2->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew2.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh2,Dummyfile, true);
    }
    
    TPZHierarquicalGrid CreateGridFrom2(GeoMesh2);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc2 = new TPZDummyFunction<STATE>(Parametricfunction2);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dt, n);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew3.txt");
        GeoMesh3->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew3.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh3,Dummyfile, true);
    }
    
    TPZHierarquicalGrid CreateGridFrom3(GeoMesh3);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc3 = new TPZDummyFunction<STATE>(Parametricfunction3);
    CreateGridFrom3.SetParametricFunction(ParFunc3);

    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh4 = CreateGridFrom3.ComputeExtrusion(t, dt, n);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew4.txt");
        GeoMesh4->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew4.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh4,Dummyfile, true);
    }
 
    return 0;
}
void Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0];
    X[1] = cos(par[0]);
    X[2] = 0.0;
}

void Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void Parametricfunction3(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = cos(par[0]);
    X[2] = par[0];
}