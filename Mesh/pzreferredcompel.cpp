/**
 * @file
 * @brief Contains the implementation of the TPZReferredCompEl methods.
 */

#include "pzreferredcompel.h"
#include "pzelctemp.h"
#include "pzintel.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"

#include "pzquad.h"
#include "TPZGeoElement.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h" 
#include "TPZGeoLinear.h"
#include "pzshapelinear.h"
#include "pzgeoquad.h"
#include "pzshapequad.h"
#include "pzgeotriangle.h"
#include "pzshapetriang.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "pzgeoprism.h"
#include "pzshapeprism.h"
#include "pzgeopyramid.h"
#include "pzshapepiram.h"
#include "pzgeotetrahedra.h"
#include "pzshapetetra.h"

#include "TPZMaterial.h"
#include "pzelmat.h"
#include "pzgeoel.h"
#include "pzcmesh.h"

#include "pzlog.h"

#include "tpzcompmeshreferred.h"

#ifdef PZDEBUG
#define DEBUG2
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompel"));
#endif

using namespace std;

template<class TVar>
void Append(TPZVec<TVar> &u1, TPZVec<TVar> &u2, TPZVec<TVar> &u12)
{
	int64_t nu1 = u1.NElements(), nu2 = u2.NElements();
	u12.Resize(nu1 + nu2);
	int64_t i;
	for (i = 0; i<nu1; i++) u12[i] = u1[i];
	for (i = 0; i<nu2; i++) u12[i + nu1] = u2[i];
}

template<class TVar>
void Append(TPZFMatrix<TVar> &u1, TPZFMatrix<TVar> &u2, TPZFMatrix<TVar> &u12)
{
	int64_t ru1 = u1.Rows(), cu1 = u1.Cols(), ru2 = u2.Rows(), cu2 = u2.Cols();
	int64_t ru12 = ru1 < ru2 ? ru2 : ru1;
	int64_t cu12 = cu1 + cu2;
	u12.Redim(ru12, cu12);
	int64_t i, j;
	for (i = 0; i<ru1; i++) for (j = 0; j<cu1; j++) u12(i, j) = u1(i, j);
	for (i = 0; i<ru2; i++) for (j = 0; j<cu2; j++) u12(i, j + cu1) = u2(i, j);
}

bool AreEqual(const TPZVec<REAL> &A, const TPZVec<REAL> &B, REAL tol) {
	if (A.NElements() != B.NElements()) return false;
	int64_t i;
	const int64_t n = A.NElements();
	for (i = 0; i < n; i++) {
		if (fabs(A[i] - B[i]) > tol) return false;
	}
	return true;
}

template< class TCOMPEL>
TPZCompEl * TPZReferredCompEl<TCOMPEL>::ReferredElement(){
	TPZCompMesh * cmesh = this->Mesh();
	TPZCompMeshReferred * refmesh = dynamic_cast<TPZCompMeshReferred*>(cmesh);
	if (!refmesh) return NULL;
	TPZCompEl * other = refmesh->ReferredEl( this->Index() );
//    if(!other) DebugStop();
	return other;
}

template< class TCOMPEL>
void TPZReferredCompEl<TCOMPEL>::Print(std::ostream & out) const{
	out << "\n" << __PRETTY_FUNCTION__ << "\n";
	TCOMPEL::Print(out);
	
	TPZCompMesh * cmesh = this->Mesh();
	TPZCompMeshReferred * refmesh = dynamic_cast<TPZCompMeshReferred*>(cmesh);
	if (refmesh){
		TPZCompEl * other = refmesh->ReferredEl( this->Index() );
		out << "My ReferredEl = " << other << "\n";
	}
	else out << "My ReferredEl = " << 0 << "\n";
	out << "end of " << __PRETTY_FUNCTION__ << "\n";
}//void

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::TPZReferredCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index):TPZRegisterClassId(&TPZReferredCompEl::ClassId),
TCOMPEL(mesh, gel,index){
	
}//method

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::TPZReferredCompEl():TPZRegisterClassId(&TPZReferredCompEl::ClassId),
TCOMPEL(){
	
}//method

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::TPZReferredCompEl(TPZCompMesh &mesh, const TPZReferredCompEl<TCOMPEL> &copy):TPZRegisterClassId(&TPZReferredCompEl::ClassId),
TCOMPEL(mesh,copy){
	
}//method

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::TPZReferredCompEl(TPZCompMesh &mesh,
											  const TPZReferredCompEl<TCOMPEL> &copy,
											  std::map<int64_t,int64_t> & gl2lcConMap,
											  std::map<int64_t,int64_t> & gl2lcElMap):
TPZRegisterClassId(&TPZReferredCompEl::ClassId),
TCOMPEL(mesh,copy,gl2lcConMap,gl2lcElMap)
{
	
}//method

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::~TPZReferredCompEl(){
	
}//method

template <>
void TPZReferredCompEl< TPZInterfaceElement >::AppendOtherSolution(TPZVec<REAL> &qsi, TPZSolVec &sol,
                                                       TPZGradSolVec &dsol, TPZFNMatrix<9,REAL> &axes)
{
    DebugStop();
}

template < class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AppendOtherSolution(TPZVec<REAL> &qsi, TPZSolVec &sol,
													   TPZGradSolVec &dsol, TPZFNMatrix<9,REAL> &axes)
{
	TCOMPEL * other = dynamic_cast<TCOMPEL *> (this->ReferredElement());
	if (!other) return;
	
	TPZSolVec ThisSol(sol);
	TPZGradSolVec ThisDSol(dsol);
	
    TPZMaterialData otherdata;
    other->InitMaterialData(otherdata);
    other->ComputeShape(qsi,otherdata);
    other->AddSolution(qsi,otherdata);
//    TPZSolVec OtherSol;
//    TPZGradSolVec OtherDSol;
    TPZGradSolVec OtherDSol2;
//    TPZFNMatrix<9> otheraxes(3,3,0.);
//    other->ComputeSolution(qsi, OtherSol, OtherDSol, otheraxes);
    int64_t numbersol = sol.size();
    OtherDSol2.resize(numbersol);
    for (int64_t is=0; is<numbersol; is++) {
        if(sol[is].NElements()){
            AdjustSolutionDerivatives(otherdata.dsol[is],otherdata.axes,OtherDSol2[is],axes);
        }
        else if(otherdata.sol[is].NElements()){
            OtherDSol2[is] = otherdata.dsol[is];
//            OtherDSol2[is] = OtherDSol[is];
            axes = otherdata.axes;
        }
        ::Append(ThisSol[is],otherdata.sol[is],sol[is]);
        ::Append(ThisDSol[is],OtherDSol2[is],dsol[is]);
    }
}

template < class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AppendOtherSolution(TPZVec<REAL> &qsi, TPZSolVec &sol)
{
    TCOMPEL * other = dynamic_cast<TCOMPEL *> (this->ReferredElement());
    if (!other) return;
    
    TPZSolVec ThisSol(sol);
    
//    TPZSolVec OtherSol;
//    TPZGradSolVec OtherDSol;
//    TPZFNMatrix<9> otheraxes(3,3,0.);
//    other->ComputeSolution(qsi, OtherSol, OtherDSol, otheraxes);
    TPZMaterialData data;
    other->ComputeSolution(qsi,data);
    int64_t numbersol = sol.size();
    for (int64_t is=0; is<numbersol; is++) {
        ::Append(ThisSol[is],data.sol[is],sol[is]);
    }
}

template < class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AppendOtherSolution(TPZVec<REAL> &qsi, TPZSolVec &sol,
													   TPZGradSolVec&dsol, const TPZFNMatrix<9,REAL> &axes)
{
    TCOMPEL * other = dynamic_cast<TCOMPEL *> (this->ReferredElement());
	if (!other) return;
	
	TPZSolVec ThisSol(sol);
	TPZGradSolVec ThisDSol(dsol);
    int64_t numbersol = sol.size();

	TPZGradSolVec OtherDSol2(numbersol);
    TPZMaterialData data;
    other->ComputeSolution(qsi,data);
    for (int64_t is=0; is<numbersol; is++) {
        if(sol[is].NElements()){
            AdjustSolutionDerivatives(data.dsol[is],data.axes,OtherDSol2[is],axes);
        }
        else if(data.sol[is].NElements()){
            OtherDSol2[is] = data.dsol[is];
            //axes = otheraxes;
        }
        ::Append(ThisSol[is],data.sol[is],sol[is]);
        ::Append(ThisDSol[is],OtherDSol2[is],dsol[is]);
    }
}

template < class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AppendOtherSolution(TPZVec<REAL> &qsi,
                                                       TPZVec<REAL> &normal,
                                                       TPZSolVec &leftsol, TPZGradSolVec &dleftsol, TPZFNMatrix<9,REAL> &leftaxes,
                                                       TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFNMatrix<9,REAL> &rightaxes){
    TCOMPEL * other = dynamic_cast<TCOMPEL *> (this->ReferredElement());
	if (!other) return;
	
	TPZSolVec ThisLeftSol(leftsol), ThisRightSol(rightsol);
	TPZGradSolVec ThisDLeftSol(dleftsol), ThisDRightSol(drightsol);
	
	TPZSolVec OtherLeftSol(0), OtherRightSol(0);
    TPZManVector<REAL,3> OtherNormal(0);
	TPZGradSolVec OtherDSol2(0), OtherDLeftSol(0), OtherDLeftSol2(0), OtherDRightSol(0), OtherDRightSol2(0);
	TPZFNMatrix<9> OtherLeftAxes(3,3,0.), OtherRightAxes(3,3,0.);
    
    TPZMaterialData dataleft, dataright;
    this->InitMaterialData(dataleft);
    this->InitMaterialData(dataright);
    
	other->ComputeSolution(qsi, OtherNormal,
						   dataleft,
						   dataright);
	
	if (OtherLeftSol.NElements() || OtherRightSol.NElements()) {//it means other has solution left/right
		if (normal.NElements() && OtherNormal.NElements()){ //then both element must have same normal
			if ( !AreEqual(normal,OtherNormal) ) {
				PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
			}
		}
		if (normal.NElements() == 0) {//however, this may be a interpolationspace and other is an interface.
			normal = OtherNormal;//Then OtherNormal is the corret value
		}//if (normal.NElements() == 0)
	}//if other has solution
	
	int64_t numbersol = ThisLeftSol.size();
    for (int64_t is=0; is<numbersol; is++) {
        if(dataleft.sol.NElements()){
            AdjustSolutionDerivatives(OtherDLeftSol[is],OtherLeftAxes,OtherDLeftSol2[is],dataleft.axes);
        }
        else if(OtherLeftSol.NElements()){
            OtherDLeftSol2[is] = OtherDLeftSol[is];
            dataleft.axes = OtherLeftAxes;
        }
	
        if(dataright.sol.NElements()){
            AdjustSolutionDerivatives(OtherDRightSol[is],OtherRightAxes,OtherDRightSol2[is],dataright.axes);
        }
        else if(OtherRightSol.NElements()){
            OtherDRightSol2[is] = OtherDRightSol[is];
            dataright.axes = OtherRightAxes;
        }
        ::Append(ThisLeftSol[is], OtherLeftSol[is], dataleft.sol[is]);
        ::Append(ThisDLeftSol[is], OtherDLeftSol[is], dataleft.dsol[is]);
        ::Append(ThisRightSol[is], OtherRightSol[is], dataright.sol[is]);
        ::Append(ThisDRightSol[is], OtherDRightSol[is], dataright.dsol[is]);
    }
}

template <  >
void TPZReferredCompEl< TPZCompElDisc >::SetCreateFunctions(TPZCompMesh *mesh){
	mesh->SetAllCreateFunctionsDiscontinuousReferred();
}

template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::SetCreateFunctions(TPZCompMesh *mesh){
	mesh->SetAllCreateFunctionsContinuousReferred();
}

/**
 * @brief Computes solution and its derivatives in local coordinate qsi
 * @param qsi master element coordinate
 * @param phi matrix containing shape functions compute in qsi point
 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
 * @param axes direction of the derivatives
 * @param sol finite element solution
 * @param dsol solution derivatives
 */
template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AddSolution(TPZVec<REAL> &qsi, TPZMaterialData &data)
{
    TCOMPEL::AddSolution(qsi,data);
    if(data.fShapeType != TPZMaterialData::EVecandShape && data.fShapeType != TPZMaterialData::EVecShape)
    {
        this->AppendOtherSolution(qsi,data.sol,data.dsol,data.axes);
    }
    else
    {
        this->AppendOtherSolution(qsi,data.sol);
    }
}


template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::ComputeSolution(TPZVec<REAL> &qsi,
                                                       TPZVec<REAL> &normal,
                                                       TPZMaterialData &dataleft,
                                                       TPZMaterialData &dataright){
    TCOMPEL::ComputeSolution(qsi, normal, dataleft, dataright);
	this->AppendOtherSolution(qsi, normal, dataleft.sol, dataleft.dsol, dataleft.axes, dataright.sol, dataright.dsol, dataright.axes);
}

template<class TCOMPEL>
int TPZReferredCompEl<TCOMPEL>::ClassId() const{
    return Hash("TPZReferredCompEl") ^ TCOMPEL::ClassId() << 1;
}

void AdjustSolutionDerivatives(TPZFMatrix<STATE> &dsolfrom, TPZFMatrix<REAL> &axesfrom,
                               TPZFMatrix<STATE> &dsolto, const TPZFMatrix<REAL> &axesto)
{
	TPZFNMatrix<9> axesinner, axesfromlocal;
	axesfrom.Transpose(&axesfromlocal);
	axesto.MultAdd(axesfromlocal,axesinner,axesinner,1.,0.);
	int64_t nderiv = dsolfrom.Rows();
	int64_t nstate = dsolfrom.Cols();
	dsolto.Resize(nderiv,nstate);
	int64_t id,jd,is;
	for(is=0; is<nstate; is++)
	{
		TPZManVector<STATE> dval(nderiv,0.);
		for(id=0; id<nderiv; id++)
		{
			for(jd=0; jd<nderiv; jd++)
			{
				dval[id] += dsolfrom(jd,is)*(STATE)axesinner(id,jd);
			}
		}
		for(id=0; id<nderiv; id++)
		{
			dsolto(id,is) = dval[id];
		}
	}
}



using namespace pzshape;
using namespace pzgeom;

template class TPZReferredCompEl< TPZInterfaceElement >;
template class TPZReferredCompEl< TPZCompElDisc >;
template class TPZReferredCompEl< TPZIntelGen<TPZShapePoint> >;
template class TPZReferredCompEl< TPZIntelGen<TPZShapeLinear> >;
template class TPZReferredCompEl< TPZIntelGen<TPZShapeQuad> >;
template class TPZReferredCompEl< TPZIntelGen<TPZShapeTriang> >;
template class TPZReferredCompEl< TPZIntelGen<TPZShapeCube> >;
template class TPZReferredCompEl< TPZIntelGen<TPZShapePrism> >;
template class TPZReferredCompEl< TPZIntelGen<TPZShapePiram> >;
template class TPZReferredCompEl< TPZIntelGen<TPZShapeTetra> >;

template class TPZRestoreClass<TPZReferredCompEl< TPZInterfaceElement >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZCompElDisc >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZIntelGen<TPZShapePoint> >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZIntelGen<TPZShapeLinear> >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZIntelGen<TPZShapeQuad> >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZIntelGen<TPZShapeTriang> >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZIntelGen<TPZShapeCube> >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZIntelGen<TPZShapePrism> >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZIntelGen<TPZShapePiram> >>;
template class TPZRestoreClass<TPZReferredCompEl< TPZIntelGen<TPZShapeTetra> >>;

TPZCompEl * CreateReferredPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
  	return new TPZReferredCompEl< TPZIntelGen<TPZShapePoint> >(mesh,gel,index);
}

TPZCompEl * CreateReferredLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZReferredCompEl< TPZIntelGen<TPZShapeLinear> >(mesh,gel,index);
}

TPZCompEl * CreateReferredQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZReferredCompEl< TPZIntelGen<TPZShapeQuad> >(mesh,gel,index);
}

TPZCompEl * CreateReferredTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZReferredCompEl< TPZIntelGen<TPZShapeTriang> >(mesh,gel,index);
}

TPZCompEl * CreateReferredCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZReferredCompEl< TPZIntelGen<TPZShapeCube> >(mesh,gel,index);
}

TPZCompEl * CreateReferredPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZReferredCompEl< TPZIntelGen<TPZShapePrism> >(mesh,gel,index);
}

TPZCompEl * CreateReferredPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZReferredCompEl< TPZIntelGen<TPZShapePiram> >(mesh,gel,index);
}

TPZCompEl * CreateReferredTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZReferredCompEl< TPZIntelGen<TPZShapeTetra> >(mesh,gel,index);
}

TPZCompEl * CreateReferredDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZReferredCompEl< TPZCompElDisc >(mesh,gel,index);
}

#include "pzelchdiv.h"

template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDiv<TPZShapeLinear>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDiv<TPZShapeTriang>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDiv<TPZShapeQuad>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDiv<TPZShapeCube>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDiv<TPZShapeTetra>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDiv<TPZShapePrism>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDiv<TPZShapePiram>>>;


template class TPZReferredCompEl<TPZCompElHDiv<TPZShapeLinear>>;
template class TPZReferredCompEl<TPZCompElHDiv<TPZShapeTriang>>;
template class TPZReferredCompEl<TPZCompElHDiv<TPZShapeQuad>>;
template class TPZReferredCompEl<TPZCompElHDiv<TPZShapeTetra>>;
template class TPZReferredCompEl<TPZCompElHDiv<TPZShapePrism>>;
template class TPZReferredCompEl<TPZCompElHDiv<TPZShapePiram>>;
template class TPZReferredCompEl<TPZCompElHDiv<TPZShapeCube>>;

#include "pzelchdivbound2.h"

template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDivBound2<TPZShapePoint>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDivBound2<TPZShapeLinear>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDivBound2<TPZShapeTriang>>>;
template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHDivBound2<TPZShapeQuad>>>;


template class TPZReferredCompEl<TPZCompElHDivBound2<TPZShapeTriang>>;
template class TPZReferredCompEl<TPZCompElHDivBound2<TPZShapePoint>>;
template class TPZReferredCompEl<TPZCompElHDivBound2<TPZShapeLinear>>;
template class TPZReferredCompEl<TPZCompElHDivBound2<TPZShapeQuad>>;

