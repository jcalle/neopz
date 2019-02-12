#include <pzgmesh.h>
#include <pzvec.h>
#include "SPZAcousticData.h"
//struct SPZAcousticTransientData{
//
//};

class TPZAcousticsSimulation {
public:
	void RunSimulation();
	SPZAcousticData fSimData;
protected:
	REAL CalculateDeltaT(std::map<int,REAL> elSizeMap, std::map<int,REAL> velocityMap);
};
//
//class TPZAcousticsFourierSimulation : public TPZAcousticsSimulation{
//protected:
//	SPZAcousticFourierData fFourierData;
//
//};
//
//class TPZAcousticsTransientSimulation : public TPZAcousticsSimulation{
//	SPZAcousticTransientData fTransientData;
//
//};