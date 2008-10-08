//
// C++ Interface: tpzintrulet3d
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZINTRULET3D_H
#define TPZINTRULET3D_H

#include "pzreal.h"
template<class T>
class TPZVec;

/**
Integration rule for tetrahedra

	@author Philippe R. B. Devloo <phil@fec.unicamp.br>
*/
class TPZIntRuleT3D{

    friend class TPZIntRuleList;
    
    short	 fNumInt;		// number of integration points for this object
    REAL	*fLocationKsi;	// location of the integration point Ksi
    REAL	*fLocationEta;	// location of the integration point Eta
    REAL	*fLocationZeta;	// location of the integration point Eta
    REAL	*fWeight;		// weight of the integration point

    TPZIntRuleT3D(int i = 2);
    ~TPZIntRuleT3D();

    public:
      enum {NUMINT_RULEST3D = 8};

      short NInt(){ return fNumInt;}	//return number of integration points

      void Loc(int i, TPZVec<REAL> &pos);			   //return location of the ith pot

      REAL W(int i);						//return weight for the ith point

};

#endif
