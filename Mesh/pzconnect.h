//$Id: pzconnect.h,v 1.7 2005-02-28 22:08:50 phil Exp $

// -*- c++ -*-
//HEADER FILE FOR CLASS NODE

#ifndef  PZCONNECTH
#define  PZCONNECTH

#include "pzfmatrix.h"
#include "pzstack.h"
#include <iostream>

using namespace std;

class TPZBndCond;
class TPZCompMesh;
class TPZBlock;
class TPZStream;



/// This class represents a set of shape functions associated with a computational element/side
/**
This class keeps track of information associated with an element/side such as order of interpolation
sequence number in the vector of blocks of equations
Objects of this class also contain the information necessary for constraints between shapefunctions
@ingroup interpolation
*/
class TPZConnect {
  /**Node block number*/
  int		fSequenceNumber;
  /**Number of element connected*/
  int		fNElConnected;
  /**
   * Interpolation order of the associated shape functions
   */
  int fOrder;

 public:
  /**Structure to reference dependency*/
  struct TPZDepend {
    int			fDepConnectIndex;
    TPZFMatrix	fDepMatrix;
    TPZDepend		*fNext;

    TPZDepend(int DepConnectIndex,TPZFMatrix &depmat,int ipos,int jpos, int isize, int jsize);

    TPZDepend(const TPZDepend &copy);
    TPZDepend(int connectindex);

    ~TPZDepend();
    TPZDepend *HasDepend(int DepConnectIndex);
    TPZDepend *RemoveDepend(TPZDepend *Ptr);
    void Write(TPZStream &buf);
    void Read(TPZStream &buf);
  };

 private:
  TPZDepend *fDependList;

 public:
  /**Constructor*/
  TPZConnect();
  /**Destructor*/
  ~TPZConnect();

  void operator=(const TPZConnect &con);

  /**Number of degrees of freedom associated with the object
     It needs the mesh object to access this information*/
  int NDof(TPZCompMesh &mesh); /*Cedric 30/09/98 -  14:05*/
  
  // Philippe In order to compute NDof, the Connect object needs to know the mesh

  /**Return the Sequence number of the connect object.
     If the sequence number == -1 this means that the node is unused*/
  int SequenceNumber() const;

  /**Set the sequence number for the global system of equations of the connect
     object. If the argument i==-1 this means that the node is out of use*/
  void SetSequenceNumber(int i) {fSequenceNumber = i;}

  /**
   * Set the order of the shapefunction associated with the connect
   */
  void SetOrder(int order) {
    fOrder = order;
  }

  /**
   * Access function to return the order associated with the connect
   */
  int Order() {
    return fOrder;
  }

  /**Print the information for the connect element. The mesh
     argument allows the object to identify the number of variables
     associated with it and the solution*/
  void Print(TPZCompMesh &mesh, ostream & out = cout);

  /**Initialize with zero fNElConnected*/
  void ResetElConnected() { fNElConnected = 0; }
  /**Increment fNElConnected*/
  void IncrementElConnected() { fNElConnected++; }
  /**Return fNElConnected*/
  int NElConnected() { return fNElConnected; }

  void AddDependency(int myindex, int dependindex,TPZFMatrix &depmat,int ipos,int jpos, int isize, int jsize);

  void RemoveDepend(int myindex, int dependindex);

  int DependencyDepth(TPZCompMesh &mesh);

  int HasDependency() { return fDependList != 0; }

  int CheckDependency(int nshape, TPZCompMesh *mesh, int nstate);


  TPZDepend *FirstDepend() { return fDependList; }

  /**Adds itself and the connects from which it depends to the list*/
  void AddToList(int myindex, TPZCompMesh &mesh, TPZStack<int> &connectlist);

  void SetDependenceOrder(int myindex, TPZCompMesh &mesh, int CurrentOrder,TPZVec<int> &connectlist,TPZVec<int> &DependenceOrder);

  void ExpandShape(int cind, TPZVec<int> &connectlist, TPZVec<int> &blocksize, TPZFMatrix &phi, TPZFMatrix &dphi);
  
  /**
  Save the element data to a stream
  */
void Write(TPZStream &buf, int withclassid);
  
  /**
  Read the element data from a stream
  */
void Read(TPZStream &buf, void *context);

};


/** 
@brief Associate a degree of freedom node with a boundary condition
such boundary condition can be dirichlet, point load or mixed boundary condition.
\deprecated
THIS CLASS IS DEPRECATED
*/
struct TPZConnectBC {

  TPZConnect *fConnect;
  TPZBndCond *fBC;

  TPZConnectBC() {
    fConnect = 0;
    fBC = 0;
  }
  TPZConnectBC(TPZConnect *nd,TPZBndCond *bc) {
    fConnect = nd;
    fBC = bc;
  }

  void Print(TPZCompMesh &mesh,ostream &out = cout);

};


#endif

