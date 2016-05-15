//
//  TRMOrchestra.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMOrchestra__
#define __PZ__TRMOrchestra__

#include <stdio.h>
#include "pzgmesh.h"
#include "tpzautopointer.h"

// Type of structural matrices
#include "TPZSkylineNSymStructMatrix.h"

#include "TRMSpaceOdissey.h"
#include "TRMSimulationData.h"
#include "TRMPrimalMultiphaseAnalysis.h"
#include "TRMMonolithicMultiphaseAnalysis.h"
#include "TRMFluxPressureAnalysis.h"
#include "TRMTransportAnalysis.h"


class TRMRawData;

class TRMOrchestra{
    
private:

    /** @brief Define the global geometry being used */
    TPZAutoPointer<TPZGeoMesh > fgmesh;

    /** @brief Define the space generator */
    TPZAutoPointer<TRMSpaceOdissey> fSpaceGenerator;
    
    /** @brief Define simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief Define the primal with global post-processing analysis */
    TPZAutoPointer<TRMPrimalMultiphaseAnalysis> fPrimalMultiphaseAnalysis;
    
    /** @brief Define the monolithic multiphase analysis */
    TPZAutoPointer<TRMMonolithicMultiphaseAnalysis> fMonolithicMultiphaseAnalysis;
    
    /** @brief Define the mixed system analysis */
    TPZAutoPointer<TRMFluxPressureAnalysis> fFluxPressureAnalysis;
    
    /** @brief Define the water transpor equation analysis */
    TPZAutoPointer<TRMTransportAnalysis> fWaterTransportAnalysis;
    
    /** @brief Define the oil transpor equation analysis */
    TPZAutoPointer<TRMTransportAnalysis> fOilTransportAnalysis;
    
    
protected:
    
    /** @brief Solve the initial conditions for pressure using a l2 projection */
    void InitializePressure();
    
public:

    /** @brief Default constructor */
    TRMOrchestra();

    /** @brief Default desconstructor */
    ~TRMOrchestra();
    
    /**
     * @defgroup Access Methods
     * @brief    Implements Access methods:
     * @{
     */
    
    /** @brief Set autopointer of the global geometry being used */
    void SetGMesh(TPZAutoPointer<TPZGeoMesh > gmesh)
    {
        fgmesh = gmesh;
    }
    /** @brief Get autopointer of the global geometry being used */
    TPZAutoPointer<TPZGeoMesh > GMesh()
    {
        return fgmesh;
    }
    
    /** @brief Set the space generator */
    void SetSpaceGenerator(TPZAutoPointer<TRMSpaceOdissey> &SpaceGenerator)
    {
        fSpaceGenerator = SpaceGenerator;
    }
    /** @brief Get the space generator */
    TPZAutoPointer<TRMSpaceOdissey> SpaceGenerator()
    {
        return fSpaceGenerator;
    }
    
    /** @brief Set autopointer of the simulation data */
    void SetSimulationData(TPZAutoPointer<TRMSimulationData > &SimulationData)
    {
        fSimulationData = SimulationData;
        fSpaceGenerator->SetSimulationData(SimulationData);
    }
    /** @brief Get autopointer of the simulation data */
    TPZAutoPointer<TRMSimulationData > SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set autopointer of the primal with global post-processing analysis */
    void SetPrimalMultiphaseAnalysis(TPZAutoPointer<TRMPrimalMultiphaseAnalysis > &PrimalMultiphaseAnalysis)
    {
        fPrimalMultiphaseAnalysis = PrimalMultiphaseAnalysis;
    }
    /** @brief Get autopointer of the primal with global post-processing analysis */
    TPZAutoPointer<TRMPrimalMultiphaseAnalysis > PrimalMultiphaseAnalysis()
    {
        return fPrimalMultiphaseAnalysis;
    }
    
    /** @brief Set autopointer of the monolithic multiphase analysis */
    void SetMonolithicMultiphaseAnalysis(TPZAutoPointer<TRMMonolithicMultiphaseAnalysis > &MonolithicMultiphaseAnalysis)
    {
        fMonolithicMultiphaseAnalysis = MonolithicMultiphaseAnalysis;
    }
    /** @brief Get autopointer of the monolithic multiphase analysis */
    TPZAutoPointer<TRMMonolithicMultiphaseAnalysis > MonolithicMultiphaseAnalysis()
    {
        return fMonolithicMultiphaseAnalysis;
    }
    
    /** @brief Set autopointer of the mixed system analysis */
    void SetFluxPressureAnalysis(TPZAutoPointer<TRMFluxPressureAnalysis > &FluxPressureAnalysis)
    {
        fFluxPressureAnalysis = FluxPressureAnalysis;
    }
    /** @brief Get autopointer of the mixed system analysis */
    TPZAutoPointer<TRMFluxPressureAnalysis > FluxPressureAnalysis()
    {
        return fFluxPressureAnalysis;
    }
    
    /** @brief Set autopointer of the water transpor equation analysis */
    void SetWaterTransportAnalysis(TPZAutoPointer<TRMTransportAnalysis > &WaterTransportAnalysis)
    {
        fWaterTransportAnalysis = WaterTransportAnalysis;
    }
    /** @brief Get autopointer of the water transpor equation analysis */
    TPZAutoPointer<TRMTransportAnalysis > WaterTransportAnalysis()
    {
        return fWaterTransportAnalysis;
    }
    
    /** @brief Set autopointer of the oil transpor equation analysis */
    void SetOilTransportAnalysis(TPZAutoPointer<TRMTransportAnalysis > &OilTransportAnalysis)
    {
        fOilTransportAnalysis = OilTransportAnalysis;
    }
    /** @brief Get autopointer of the oil transpor equation analysis */
    TPZAutoPointer<TRMTransportAnalysis > OilTransportAnalysis()
    {
        return fOilTransportAnalysis;
    }
    
    // @}
    
    /**
     * @defgroup Methods for run different analysis
     * @brief    Analysis creation and excecution:
     * @{
     */
    
    /** @brief Create computational meshes using space odissey */
    void CreateCompMeshes();
    
    /** @brief Create a primal analysis using space odissey */
    void CreateAnalysisPrimal();

    /** @brief Create a dual analysis using space odissey */
    void CreateAnalysisDual();
    
    /** @brief Create a dual analysis on box geometry using space odissey */
    void CreateAnalysisDualonBox();
    
    /** @brief Create a monolithic dual analysis on box geometry using space odissey */
    void CreateMonolithicAnalysis();
    
    /** @brief Create a monolithic dual analysis on box geometry using space odissey */
    void OneStepMonolithicAnalysis();
    
    /** @brief Create a monolithic dual analysis on box geometry using space odissey */
    void PostProMonolithicAnalysis();
    
    /** @brief Run the time steps set in the simulation data */
    void RunSimulation();
    
    /** @brief Run a single time step */
    void ExecuteOneTimeStep();
    
    /** @brief Computes the post processed results */
    void PostProcess();
    
    /** @brief Project an exact solution */
    void ProjectExactSolution();
    
    /** @brief Compute the production rate of the reservoir */
    void ComputeProductionRate(std::map<REAL,STATE> &RatebyPosition, STATE &Total);
    
    /** @brief exact pressure */
    static void ExactPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &pressure);
    
    /** @brief exact flux */
    static void ExactFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &flux);
    
    /** @brief exact laplacian */
    static void ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &pressure);
    
    /** @brief Compute the system of equations using transfer matrixces */
    TPZFMatrix<STATE> IntegrateResidue(TPZAutoPointer<TPZCompMesh> cmesh_multiphysics, TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer< TPZCompMesh> cmesh_pressure, TPZAutoPointer<TRMBuildTransfers> transfer);
    
    /** @brief Compute gradient of the system of equations using transfer matrixces */
    void IntegrateGradientOfResidue(TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer<TRMBuildTransfers> transfer);
    
    // @}    
    
};

#endif /* defined(__PZ__TRMOrchestra__) */
