/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticTrig methods. 
 */
#include "tpzquadratictrig.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzgeotriangle.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "tpztriangle.h"
#include "pznoderep.h"
#include "pzshapetriang.h"
#include "tpzgeoelmapped.h"

#include <iostream>
#include <iostream>
#include <cmath>

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"

#include "tpzgeomid.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

/**
 / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
 / LabMeC - FEC - UNICAMP
 / 2007
 */

template<class T>
void TPZQuadraticTrig::Shape(TPZVec<T> &param,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi)
{
	T qsi = param[0], eta = param[1];
    
	phi(0,0) = (qsi+eta-1.)*(2.*qsi+2.*eta-1.);
	phi(1,0) = qsi*(2.*qsi-1.);
	phi(2,0) = eta*(2.*eta-1.);
	phi(3,0) = -4.*qsi*(qsi+eta-1.);
	phi(4,0) =  4.*qsi*eta;
	phi(5,0) = -4.*eta*(qsi+eta-1.);
	dphi(0,0) = 1.-4.*(1.-qsi-eta);
    dphi(0,1) = -1.+4.*qsi;
    dphi(0,2) = 0.;
    dphi(0,3) = 4.*(1.-qsi-eta)-4.*qsi;
    dphi(0,4) = 4.*eta;
    dphi(0,5) = -4.*eta;
	dphi(1,0) = 1.-4.*(1.-qsi-eta);
    dphi(1,1) = 0.;
    dphi(1,2) = -1.+4.*eta;
    dphi(1,3) = -4.*qsi;
    dphi(1,4) = 4.*qsi;
    dphi(1,5) = 4.*(1.-qsi-eta)-4.*eta;
}

template<class T>
void TPZQuadraticTrig::X(TPZFMatrix<REAL> &coord, TPZVec<T>& par, TPZVec< T >& result)
{
	TPZManVector<T,3> parMap(2);
	REAL spacephi[6],spacedphi[12];
	TPZFNMatrix<6,T> phi(6,1);
	TPZFNMatrix<12,T> dphi(2,6);
	Shape(par,phi,dphi);
	for(int i = 0; i < 3; i++)
	{
		result[i] = 0.0;
		for(int j = 0; j < 6; j++) result[i] += phi(j,0)*coord(i,j);
	}
}

void TPZQuadraticTrig::Jacobian(TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv)
{
	jacobian.Resize(2,2); axes.Resize(2,3); jacinv.Resize(2,2);
	REAL spacephi[6],spacedphi[12];
	TPZFMatrix<REAL> phi(6,1,spacephi,6);
	TPZFMatrix<REAL> dphi(2,6,spacedphi,12);
	Shape(par,phi,dphi);
	jacobian.Zero();
	
	TPZFMatrix<REAL> VecMatrix(3,2,0.);
	for(int i = 0; i < 6; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			VecMatrix(j,0) += coord(j,i)*dphi(0,i);
			VecMatrix(j,1) += coord(j,i)*dphi(1,i);
		}
	}
	TPZFMatrix<REAL> axest;
	VecMatrix.GramSchmidt(axest,jacobian);
	axest.Transpose(&axes);
	
	detjac = jacobian(0,0)*jacobian(1,1)-jacobian(1,0)*jacobian(0,1);
	jacinv(0,0) =  jacobian(1,1)/detjac;
	jacinv(1,1) =  jacobian(0,0)/detjac;
	jacinv(0,1) = -jacobian(0,1)/detjac;
	jacinv(1,0) = -jacobian(1,0)/detjac;
}

TPZGeoEl *TPZQuadraticTrig::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc)
{
	if(side==6)
	{
		TPZManVector<long> nodes(3); int i;
		for (i=0;i<3;i++) nodes[i] = orig->SideNodeIndex(side,i);
		long index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoBlendElement(ETriangle,nodes,bc,index);
		int iside;
		for (iside = 0; iside <6; iside++)
		{
			TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::ContainedSideLocId(side,iside)));
		}
		TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	
	else if(side>-1 && side<3)
	{
		TPZManVector<long> nodeindexes(1);
		nodeindexes[0] = orig->SideNodeIndex(side,0);
		long index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoBlendElement(EPoint,nodeindexes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	
	else if(side > 2 && side < 6)
	{
		TPZManVector<long> nodes(2);
		nodes[0] = orig->SideNodeIndex(side,0);
		nodes[1] = orig->SideNodeIndex(side,1);
		long index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoBlendElement(EOned,nodes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::ContainedSideLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::ContainedSideLocId(side,1)));
		TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
        gel->BuildBlendConnectivity();
		return gel;
	}
	
	else PZError << "TPZGeoTriangle::CreateBCGeoEl has no bc.\n";
	return 0;
}


/**
 * Creates a geometric element according to the type of the father element
 */
TPZGeoEl *TPZQuadraticTrig::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											 TPZVec<long>& nodeindexes,
											 int matid,
											 long& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

/// create an example element based on the topology
/* @param gmesh mesh in which the element should be inserted
 @param matid material id of the element
 @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
 @param size (in) size of space where the element should be created
 */
#include "tpzchangeel.h"

void TPZQuadraticTrig::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    TPZManVector<REAL,3> co(3),shift(3),scale(3);
    TPZManVector<long,3> nodeindexes(3);
    for (int i=0; i<3; i++) {
        scale[i] = size[i]/3.;
        shift[i] = 1./2.+lowercorner[i];
    }
    
    for (int i=0; i<NCornerNodes; i++) {
        ParametricDomainNodeCoord(i, co);
        for (int j=0; j<co.size(); j++) {
            co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
        }
        nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
    }
    long index;
    CreateGeoElement(gmesh, ETriangle, nodeindexes, matid, index);
    TPZGeoEl *gel = gmesh.Element(index);
    int nsides = gel->NSides();
    for (int is=0; is<nsides; is++) {
        gel->SetSideDefined(is);
    }
    gel = TPZChangeEl::ChangeToQuadratic(&gmesh, index);
    for (int node = gel->NCornerNodes(); node < gel->NNodes(); node++) {
        TPZManVector<REAL,3> co(3);
        gel->NodePtr(node)->GetCoordinates(co);
        for (int i=0; i<3; i++) {
            co[i] += (0.2*rand())/RAND_MAX - 0.1;
        }
        gel->NodePtr(node)->SetCoord(co);
    }
}


//void TPZQuadraticTrig::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
//{
//    if(node > this->NNodes)
//    {
//        DebugStop();
//    }
//    nodeCoord.Resize(Dimension, 0.);
//    switch (node) {
//        case (0):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 0.;
//            break;
//        }
//        case (1):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 0.;
//            break;
//        }
//        case (2):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 1.;
//            break;
//        }
//        case (3):
//        {
//            nodeCoord[0] = 0.5;
//            nodeCoord[1] = 0.0;
//            break;
//        }
//        case (4):
//        {
//            nodeCoord[0] = 0.5;
//            nodeCoord[1] = 0.5;
//            break;
//        }
//        case (5):
//        {
//            nodeCoord[0] = 0.0;
//            nodeCoord[1] = 0.5;
//            break;
//        }
//        default:
//        {
//            DebugStop();
//            break;
//        }
//    }
//}

///CreateGeoElement -> TPZQuadraticTrig

template<>
int TPZGeoElRefPattern<TPZQuadraticTrig>::ClassId() const {
	return TPZGEOELEMENTQUADRATICTRIANGLEID;
}

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticTrig>, TPZGEOELEMENTQUADRATICTRIANGLEID>;

template class TPZGeoElRefLess<TPZQuadraticTrig>;
template class pzgeom::TPZNodeRep<6,TPZQuadraticTrig>;

