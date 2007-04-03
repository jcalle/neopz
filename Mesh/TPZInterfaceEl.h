// -*- c++ -*-

//$Id: TPZInterfaceEl.h,v 1.41 2007-04-03 12:30:26 tiago Exp $

#ifndef ELEMINTERFACEHH
#define ELEMINTERFACEHH

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
class TPZCompElDisc;

/// This class computes the contribution over an interface between two discontinuous elements
/**
@ingroup CompElement
*/
class TPZInterfaceElement : public TPZCompEl {

   private :

  /**
   * Element the left of the normal a interface
   */
  TPZCompElSide fLeftElSide;

  /**
   * element the right of the normal a interface
   */
  TPZCompElSide fRightElSide;

  /**
   * Normal to the face element
   */
  TPZManVector<REAL,3> fNormal;

  /**
   * Keep track of the connects of the element
   */
  //TPZVec<TPZConnect *> fConnectL, fConnectR;

 /**
   * Keep track of the connect indexes
   */
//  int fConnectIndexL[NL], fConnectIndexR[NR];

  /**
   * Geometric element to which this element refers
   */
//  TPZGeoEl *fReference;

  /**
   * Material object of this element
   */
  int fMaterialId;//this variable can be gotten of the element of associated volume

  void DecreaseElConnected();
  void IncrementElConnected();

 /**
  * Extract connects from element el.
  */
  void GetConnects(TPZCompElSide &elside, TPZVec<TPZConnect*> &connects, TPZVec<int> &connectindex);

  /**
   * Compute shape functions to an interpolated element. Used in case one neighbour is TPZInterpolatedElement.
   */
  void ComputeShape(TPZInterpolatedElement* intel, TPZFMatrix &phix, TPZFMatrix &dphix,
                    TPZFMatrix &axes, TPZVec<REAL> &IntPoint );

 public:

  enum CalcStiffOptions{ENone = -1, EStandard /*Deprecated*/ = 0, EPenalty, EContDisc,EReferred};

  /**
   * For CloneInterface usage. Normal is not recomputed, but copied.
   * Only for disconitnuous neighbours.
   */
  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompEl *left,TPZCompEl *right, const TPZVec<REAL> & normal);

  /**
   * Construtor para o elemento descontinuo.
   */
//  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompEl *left,TPZCompEl *right);

  /** Constuctor to continuous and/or discontinuous neighbours.
   */
  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompEl *left,TPZCompEl *right, int leftside, int rightside);

  /**
   * Copy constructor.
   */
  TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy);


  /**
   * Clone constructor to a patch mesh
   * @param mesh reference to a clone mesh
   * @param copy element to be copied
   * @param gl2lcIdx map between global(original) and local (patch) connect indexes
   */
  TPZInterfaceElement(TPZCompMesh &mesh,
                      const TPZInterfaceElement &copy,
                      std::map<int,int> &gl2lcConIdx,
                      std::map<int,int> &gl2lcElIdx);


  /**
   *
   */
  TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy, int &index);

  /**
   * Empty constructor.
   */
  TPZInterfaceElement();

  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index);
//   TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy, TPZVec<int> &destindex,int &index);

  ~TPZInterfaceElement();

  void SetLeftRightElements(TPZCompElSide & left, TPZCompElSide & right);

  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    return new TPZInterfaceElement(mesh, *this);
  }

  /**
   * @see class TPZCompEl
   */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> &gl2lcConMap, std::map<int,int> &gl2lcElMap) const
  {
    return new TPZInterfaceElement(mesh, *this, gl2lcConMap,gl2lcElMap);
  }


  TPZCompEl * CloneInterface(TPZCompMesh &aggmesh,int &index, TPZCompElDisc * left, TPZCompElDisc * right) const;

//  TPZAutoPointer<TPZMaterial> Material() const
//  {
//    return TPZAutoPointer<TPZMaterial> (fMaterial);
//  }

//  void SetMaterial(TPZAutoPointer<TPZMaterial> mat) { fMaterial = mat;}

  /**
   * it identifies the elements of left and right volume of the interface
   */
  void VolumeEls(TPZCompEl &thirdel);

  void GetTransformsLeftAndRight(TPZTransform &tl,TPZTransform &tr);

  /**
   * it returns the right element from the element interface
   */
  TPZCompEl *RightElement() {return fRightElSide.Element();}

  /**
   * it returns the left element from the element interface
   */
  TPZCompEl *LeftElement() {return fLeftElSide.Element();}

  /**
   * it returns the normal one to the face from the element
   */
  void Normal(TPZVec<REAL> &normal);

/*   void SetNormal(TPZVec<REAL> &normal); */

  /**
   * it returns the number from connectivities of the element
   */
  int NConnects();

  /**
   * it returns the number from connectivities of the element related to right neighbour
   */
  int NRightConnects();

  /**
   * it returns the number from connectivities of the element related to left neighbour
   */
  int NLeftConnects();

  /**
   * Its return the connects of the left and right element associates
   */
  int ConnectIndex(int i);

  /**
   * This function should not be called
   */
  void SetConnectIndex(int node, int index);

  /**
   * it returns the dimension from the element interface
   */
  int Dimension() const {
     return this->Reference()->Dimension();
  }

  /**
   * Type of the element
   */
  MElementType Type() { return EInterface; }

  /**
   * it returns the shapes number of the element
   * the associated space of interpolation is gotten
   * of the elements left and right
   */
  int  NShapeF() {return 0;}

  /**
   * Loads the solution within the internal data structure of the element
   * Is used to initialize the solution of connect objects with dependency
   * Is also used to load the solution within SuperElements*/
  virtual void LoadSolution(){
    //NOTHING TO BE DONE HERE
  }

  /**
   * CalcStiff computes the element stiffness matrix and right hand side
   * @param ek element matrix
   * @param ef element right hand side
   */
  virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * CalcResidual only computes the element residual
   * @param ef element residual
   */
  virtual void CalcResidual(TPZElementMatrix &ef);

  /**
   * Standard CalcStiff of interface element
   * @param ek element matrix
   * @param ef element right hand side
   */
  void CalcStiffStandard(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * Standard CalcResidual of interface element
   * @param ek element matrix
   * @param ef element right hand side
   */
  void CalcResidualStandard(TPZElementMatrix &ef);

  /**
   * CalcStiff with penalty term based on left and right p order and
   * interface inner radius.
   * @param ek element matrix
   * @param ef element right hand side
   * @author Paulo B�ing & Igor Mozolevski
   * @since April 01, 2004
   */
  void CalcStiffPenalty(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * CalcStiff for meshes who combine continuous and discontinuous
   * elements. It was not necessary to separate this implementation
   * from Standard and Penalty implementations.
   * @param ek element matrix
   * @param ef element right hand side
   * @since March 01, 2005
   */
  void CalcStiffContDisc(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * CalcStiff for meshes who combine continuous and discontinuous
   * elements and use referred meshes.
   * It was not necessary to separate this implementation
   * from Standard and Penalty implementations.
   * @param ek element matrix
   * @param ef element right hand side
   * @since March 01, 2005
   */
  void CalcStiffReferred(TPZElementMatrix &ek, TPZElementMatrix &ef);
  /**
   * gCalcStiff = 1 means standard CalcStiff
   * gCalcStiff = 2 means CalcStiff with penalty
   * gCalcStiff = 3 means the mesh has continuous and discontinuous elements combined
   * gCalcStiff = 4 means the mesh has referred meshes
   */
  static int gCalcStiff;

  static void SetCalcStiffStandard(){ TPZInterfaceElement::gCalcStiff = EStandard; }

  static void SetCalcStiffPenalty(){ TPZInterfaceElement::gCalcStiff = EPenalty; }

  static void SetCalcStiffContDisc(){ TPZInterfaceElement::gCalcStiff = EContDisc; }

  static void SetCalcStiffReferred(){ TPZInterfaceElement::gCalcStiff = EReferred; }

 /**
  * Computes solution and its derivatives in the local coordinate qsi.
  * @param [in] qsi master element coordinate
   * @param [out] leftsol left finite element solution
   * @param [out] rightsol right finite element solution
   * @param [out] dleftsol left solution derivatives
   * @param [out] drightsol right solution derivatives
   * @param [out] leftaxes axes associated with the derivative of the left element
   * @param [out] rightaxes axes associated with the derivative of the right element
  */
virtual void ComputeSolution(TPZVec<REAL> &qsi,
                               TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                               TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes);

 /**
   * Computes solution and its derivatives in the local coordinate qsi.
   * @param [in] qsi master element coordinate
   * @param [out] sol finite element solution at the interface element
   * @param [out] leftsol left finite element solution
   * @param [out] rightsol right finite element solution
   * @param [out] dsol solution derivatives at the interface element
   * @param [out] dleftsol left solution derivatives
   * @param [out] drightsol right solution derivatives
   * @param [out] axes axes associated with the derivative of the solution of the interface element
   * @param [out] leftaxes axes associated with the derivative of the left element
   * @param [out] rightaxes axes associated with the derivative of the right element
  */
virtual void ComputeSolution(TPZVec<REAL> &qsi,
                             TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes,
                             TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                             TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes)
{
  // usually interface elements have no associated solution
  return ComputeSolution(qsi,leftsol,dleftsol,leftaxes,rightsol,drightsol,rightaxes);
}

  /**
  * Computes solution and its derivatives in local coordinate qsi
  * @param qsi master element coordinate
  * @param phi matrix containing shape functions compute in qsi point
  * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
  * @param [in] axes indicating the direction of the derivatives
  * @param sol finite element solution
  * @param dsol solution derivatives
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                               TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol){
    //NOTHING TO BE DONE HERE - Interface elements have no solution associated
  }

 /**
  * Computes solution and its derivatives in the local coordinate qsi.
  * @param qsi master element coordinate
  * @param sol finite element solution
  * @param dsol solution derivatives
  * @param axes axes associated with the derivative of the solution
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi, 
                               TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes){
    //NOTHING TO BE DONE HERE - Interface elements have no solution associated
  }

  void VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet);

  /**
   * Print attributes of the object
   */
  void Print(std::ostream &out = std::cout);

  /**
   * it verifies the existence of interfaces associates
   * with the side of an element
   * case to interface should exist and exists only a returns 1
   * case to interface should not exist and does not exist returns 1
   * otherwise returns 0
   */
  static int ExistInterfaces(TPZCompElSide &comp);

  //it returns the normal one to the face from element
  void NormalToFace(TPZVec<REAL> &normal/*,int leftside*/);

  static int FreeInterface(TPZCompMesh &cmesh);

  /**
   * reproduz na malha aglomerada aggmesh uma copia da interface da malha fina
   */
  void CloneInterface(TPZCompMesh *aggmesh);

  static int main(TPZCompMesh &cmesh);

  void EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
		     TPZVec<REAL> &errors, TPZBlock * /*flux */);

  /**
   * ComputeError computes the element error estimator
  */
  virtual void ComputeError(int errorid, TPZVec<REAL> &errorL, TPZVec<REAL> &errorR);

  /**
   * Integrate a variable over the element.
   */
   virtual void Integrate(int variable, TPZVec<REAL> & value);

  void EvaluateInterfaceJumps(TPZVec<REAL> &errors);

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
  virtual int ClassId() const;
  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);

};

//Acessar com -> TPZGeoElXXd::SetCreateFunction(createInterfaceXXEl);
#endif

