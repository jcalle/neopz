/* class that defines the default refinement of the prism element */
// -*- c++ -*-

#ifndef TPZREFPRISMH
#define TPZREFPRISMH

#include "pzstack.h"
class TPZGeoEl;
class TPZGeoElSide;
class TPZTransform;

namespace pzrefine {

/// implements the uniform refinement of a geometric prism element
class TPZRefPrism {

public:

	enum{NSubEl = 8};

	static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
	static void MidSideNodeIndex(TPZGeoEl *gel,int side,int &index);
	static void NewMidSideNode(TPZGeoEl *gel,int side,int &index);
	static void GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
	static int NSideSubElements(int side);
	//static int NSideSubElements(int side);
	static TPZTransform GetTransform(int side,int son);
	static int FatherSide(int side,int son);
	//static int NSubElements();
};

};
#endif
