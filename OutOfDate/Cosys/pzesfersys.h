//HEADER FILE FOR CLASS COSYS

#ifndef ESFERSYSTHPP
#define ESFERSYSTHPP

#include "pzcosys.h"
#include "pzreal.h"

/**
 * @ingroup cosys
 * @brief Defines the esferical coordinate system
 */
class TPZEsfersys : public TPZCosys {

 public:

  /**
   * Default empty constructor 
   */
  TPZEsfersys();

  /**
   * Creates a new object from one reference object and new origin 
   * @param num index of the coordinate system
   * @param ref reference coordinate system
   */
  TPZEsfersys(int num,  TPZCartsys* ref = NULL);

  /**
   * Destructor 
   */
  ~TPZEsfersys() {;}

  /**
   * Transform point coordinates form current system to the reference system 
   * @param point point to transfer to reference coordinate system     
   */
  void ToReference(TPZVec<REAL> &point);

  /**
   * Transform the point coordinates from the reference system to the current system 
   * @param point  point to transfer from reference coordinate system
   */
  void FromReference(TPZVec<REAL> &point);

  /**
   * Return the coordinate system type 
   */
  int Type() { return esferic;}

  /**
   * Verify if the difference between two nodes are greater than PI 
   * @param points points to verify the angle
   */
  void VerifyRange(TPZFMatrix &point);

  /**
   * Calculates the transformation gradient Gradx given Gradient of  the point X 
   * @param X point that want to transform gradient
   * @param GradX gradient of point X
   * @param x will receive point X in terms of the destintation coordinate system
   * @param Gradx will receive the tranformed gradient
   * @param dest coordinate system to transform the gradient   
   */
  void TransformGradient(TPZVec<REAL> &X, TPZFMatrix &GradX, TPZVec<REAL> &x, TPZFMatrix &Gradx, TPZCosys *dest = 0);
};

#endif
