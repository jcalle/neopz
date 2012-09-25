/*
 *  TPZMohrCoulomb.h
 *  FEMPZ
 *
 *  Created by Diogo Cecilio on 5/4/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


/* Generated by Together */// $Id: TPZMohrCoulomb.h,v 1.2 2010-06-11 22:12:14 diogo Exp $

#ifndef TPZMOHRCOULOMB_H
#define TPZMOHRCOULOMB_H

#include "pzlog.h"
#include "TPZPlasticStep.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzvec_extras.h"
#include "pzsave.h"
#include "TPZPlasticStepID.h"

#ifdef LOG4CXX_PLASTICITY
static LoggerPtr loggerMohrCoulomb(Logger::getLogger("MCC"));
#endif

#define MOHRCOULOMBPARENT TPZPlasticStep<TPZYCMohrCoulomb, TPZThermoForceA, TPZElasticResponse>


class TPZMohrCoulomb : public MOHRCOULOMBPARENT, public TPZSaveable  {
	
public:
	
	enum {NYield = TPZYCMohrCoulomb::NYield};
	
public:
	
    TPZMohrCoulomb():MOHRCOULOMBPARENT()
    {
		fMaterialTensionSign  = 1; // internally in this material tension is negative
		fInterfaceTensionSign =  1; // by default
    }
	
    TPZMohrCoulomb(const TPZMohrCoulomb & source):MOHRCOULOMBPARENT(source)
    {
    }
	
    TPZMohrCoulomb & operator=(const TPZMohrCoulomb & source)
    {
		MOHRCOULOMBPARENT::operator=(source);		
		return *this;
    }
	
	virtual const char * Name() const
	{
		return "TPZMohrCoulomb";	
	}
	
    static void ConventionalConcrete(TPZMohrCoulomb & material)
	{
    	REAL pi = M_PI;
	    REAL cohesion = 11.2033; //yield- coesao inicial correspondeno a fck igual 32 Mpa
	    REAL phi =  20./180. * pi; //phi=20
	    REAL hardening = 1000.; //Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao
	    REAL young = 20000.;
	    REAL poisson = 0.2;
        material.fYC.SetUp(phi);
		material.fTFA.SetUp(cohesion, hardening);
		material.fER.SetUp(young, poisson);
	}
    
    static void TaludeMaterial(TPZMohrCoulomb & material)
	{
    	REAL pi = M_PI;
        REAL cohesion = 50.; //yield- coesao inicialem KPa
	    REAL phi =  20./180. * pi; //phi=20
	    REAL hardening = 10.; //Modulo de hardening da coesao equivante 0.01 Mpa a cada 0.1% de deformacao
	    REAL young = 20000.;//E em KPa
	    REAL poisson = 0.49;
        material.fYC.SetUp(phi);
		material.fTFA.SetUp(cohesion, hardening);
		material.fER.SetUp(young, poisson);
	}
	
	void SetUp(REAL & cohesion, REAL & phi, REAL & hardening, REAL &young, REAL &poisson)
	{
		MOHRCOULOMBPARENT::fYC.SetUp(phi);
		MOHRCOULOMBPARENT::fTFA.SetUp(cohesion, hardening);
		MOHRCOULOMBPARENT::fER.SetUp(young, poisson);
	}
	
	virtual void Print(std::ostream & out) const
	{
		out << "\n" << this->Name();
		out << "\n Base Class Data:\n";
		MOHRCOULOMBPARENT::Print(out);		
	}
	
	virtual int ClassId() const
	{	
		return TPZMOHRCOULOMB_ID;	
	}
	
	virtual void Write(TPZStream &buf, int withclassid)
	{
		TPZSaveable::Write(buf, withclassid);
				
		buf. Write(&fYC.fPhi, 1);
		
        buf. Write(&fER.fLambda, 1);
        buf. Write(&fER.fMu, 1);	
        
        buf.Write(&fTFA.fSigmaYield0, 1);
        buf.Write(&fTFA.fK, 1);

        buf. Write(&fResTol, 1);
        buf. Write(&fIntegrTol, 1);
        buf. Write(&fMaxNewton, 1);
        buf. Write(&fMinLambda, 1);
		
        buf. Write(&fN.fEpsT.fData[0], 6);
        buf. Write(&fN.fEpsP.fData[0], 6);
        buf. Write(&fN.fAlpha, 1);
	
	   // fPlasticMem does not need to be stored
		
	}

	virtual void Read(TPZStream &buf, void *context)
	{
	   TPZSaveable::Read(buf, context);

        buf. Read(&fYC.fPhi, 1);
		
        buf. Read(&fER.fLambda, 1);
        buf. Read(&fER.fMu, 1);	
		
        buf. Read(&fTFA.fSigmaYield0, 1);
        buf. Read(&fTFA.fK, 1);	
		
        buf. Read(&fResTol, 1);
        buf. Read(&fIntegrTol, 1);
        buf. Read(&fMaxNewton, 1);
        buf. Read(&fMinLambda, 1);
		
        buf. Read(&fN.fEpsT.fData[0], 6);
        buf. Read(&fN.fEpsP.fData[0], 6);
        buf. Read(&fN.fAlpha, 1);
		
        fPlasticMem.Resize(0);
	}	
	
		
public:

};


#endif //TPZMohrCoulomb_H
