//
//  TRMPhaseTransport.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPhaseTransport.h"

TRMPhaseTransport::TRMPhaseTransport() : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>()
{
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
    /** @brief material dimension */
    fdimension = 0;
    
}

TRMPhaseTransport::TRMPhaseTransport(int matid, int dimension) : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>(matid)
{
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
    /** @brief material dimension */
    fdimension = dimension;
    
}


TRMPhaseTransport::TRMPhaseTransport(const TRMPhaseTransport &mat) : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>(mat)
{
    this->fSimulationData   = mat.fSimulationData;
    this->fTransfer         = mat.fTransfer;
    this->fdimension        = mat.fdimension;
}

TRMPhaseTransport::~TRMPhaseTransport()
{
    
}

void TRMPhaseTransport::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TRMPhaseTransport::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TRMPhaseTransport::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMPhaseTransport::VariableIndex(const std::string &name) {
    if (!strcmp("sw", name.c_str())) return 0;
    if (!strcmp("so", name.c_str())) return 1;
    if (!strcmp("sg", name.c_str())) return 2;
    if (!strcmp("id", name.c_str())) return 3;

    return TPZMatWithMem::VariableIndex(name);
}

int TRMPhaseTransport::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 1; // Vector
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
    }
    return TPZMatWithMem::NSolutionVariables(var);
}

void TRMPhaseTransport::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            Solution_ab(datavec, var, Solout);
        }
            break;
        case 3:
        {
            Solution_abc(datavec, var, Solout);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}

// Jacobian contribution
void TRMPhaseTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            Contribute_ab(datavec, weight, ek, ef);
        }
            break;
        case 3:
        {
            Contribute_abc(datavec, weight, ek, ef);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}


void TRMPhaseTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            Contribute_ab(datavec, weight, ef);
        }
            break;
        case 3:
        {
            Contribute_abc(datavec, weight, ef);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}


// Two phase case

void TRMPhaseTransport::Solution_ab(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int sb_a    = 0;
    REAL sa = datavec[sb_a].sol[0][0];
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = sa;
        }
            break;
        case 1:
        {
            Solout[0] = 1.0-sa;
        }
            break;
        case 3:
        {
            Solout[0] = this->Id();
        }
            break;
        default:
        {
            DebugStop();
        }
    }
    
}

void TRMPhaseTransport::Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    
    int sb_a    = 0;
    
    TPZFNMatrix<100,STATE> phi_ss       = datavec[sb_a].phi;
    REAL sa = datavec[sb_a].sol[0][0];
    
    int nphis_a     = phi_ss.Rows();
    int firsts_a    = 0;
    
    // Time
    STATE dt = fSimulationData->dt();
    
    // Get the pressure at the integrations points
    long global_point_index = datavec[sb_a].intGlobPtIndex;
    TRMPhaseMemory & memory = GetMemory()[global_point_index];

    REAL sw      = memory.sa();
    REAL sw_n    = memory.sa_n();
    
    REAL p      = memory.p_avg();
    REAL p_n    = memory.p_avg_n();
    
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    
    // Rock parameters form point memory
    REAL phi_0;
    K       = memory.K_0();
    Kinv    = memory.Kinv_0();
    phi_0   = memory.phi_0();
    
    //  Computing closure relationship at given average values
    TPZManVector<STATE, 10> v(nvars), v_n(nvars);
    v[0]        = p;
    v[1]        = sw;
    v_n[0]      = p_n;
    v_n[1]      = sw_n;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho_w,rho_o,Bw_n,Bw;
    fSimulationData->AlphaProp()->Density(rho_w, v);
    fSimulationData->AlphaProp()->B(Bw, v);
    fSimulationData->AlphaProp()->B(Bw_n, v_n);
    
    REAL l_dr   = memory.lambda();
    REAL mu_dr  = memory.mu();
    REAL alpha  = memory.alpha();
    REAL Se = memory.S_e();
    
    REAL phi = phi_0;
    REAL phi_n = phi_0;
//    REAL Ss = 0.0;
    
    for (int is = 0; is < nphis_a; is++)
    {
        
        ef(is + firsts_a) += weight * (1.0/dt) * (phi_n*(sw_n/Bw_n[0]) - phi*(sw/Bw[0]) )* phi_ss(is,0);
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(is + firsts_a, js + firsts_a) += weight * (1.0/dt) * phi_n*(1.0/Bw_n[0]) * phi_ss(js,0) * phi_ss(is,0);
        }
        
    }
    
}

void TRMPhaseTransport::Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute_ab(datavec, weight, ek_fake, ef);
    return;
}


// Three phase case

void TRMPhaseTransport::Solution_abc(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int sb_a    = 0;
    int sb_b    = 1;
    REAL sa = datavec[sb_a].sol[0][0];
    REAL sb = datavec[sb_b].sol[0][0];
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = sa;
        }
            break;
        case 1:
        {
            Solout[0] = sb;
        }
            break;
        case 2:
        {
            Solout[0] = 1.0 - sa - sb;
        }
            break;
        default:
        {
            DebugStop();
        }
    }
    
}

void TRMPhaseTransport::Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int nvars = 4; // {p,sa,sb,t}
    
    int sb_a    = 0;
    int sb_b    = 1;
    
    TPZFNMatrix<100,STATE> phi_ssa       = datavec[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb       = datavec[sb_b].phi;
    REAL sa = datavec[sb_a].sol[0][0];
    REAL sb = datavec[sb_b].sol[0][0];
    
    int nphis_a     = phi_ssa.Rows();
    int nphis_b     = phi_ssb.Rows();
    int firsts_a    = 0;
    int firsts_b    = firsts_a + nphis_a;
    
    // Time
    STATE dt = fSimulationData->dt();
    
    // Get the pressure at the integrations points
    long global_point_index = datavec[sb_a].intGlobPtIndex;
    TRMPhaseMemory &point_memory = GetMemory()[global_point_index];
    REAL p_avg_n    = point_memory.p_avg_n();
    REAL p_avg      = point_memory.p_avg();
    
//    REAL sa_avg_n    = point_memory.sa_n();
//    REAL sa_avg      = point_memory.sa();
//
//    REAL sb_avg_n    = point_memory.sb_n();
//    REAL sb_avg      = point_memory.sb();
    
    
    
    //  Average values
    
    //  Computing closure relationship at given average values
    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_avg_n;
    v[1] = sa;
    v[2] = sb;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho_a,rho_b,l;
    fSimulationData->AlphaProp()->Density(rho_a, v);
    fSimulationData->BetaProp()->Density(rho_b, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi(nvars,0.0);
    //    fSimulationData->Map()->Kappa(datavec[ub].x, K, Kinv, v);
    //    fSimulationData->Map()->phi(datavec[ub].x, phi, v);
    
    // Rock parameters form point memory
    REAL phi_0;
    K       = point_memory.K_0();
    Kinv    = point_memory.Kinv_0();
    phi_0   = point_memory.phi_0();
    phi[0] = phi_0;
    
    
    // Integration point contribution
    
    if(! fSimulationData->IsCurrentStateQ()){
        v[0] = p_avg;
        v[1] = sa;
        v[2] = sb;
        fSimulationData->AlphaProp()->Density(rho_a, v);
        fSimulationData->BetaProp()->Density(rho_b, v);
        fSimulationData->Map()->phi(datavec[sb_a].x, phi, v);
        
        
        for (int is = 0; is < nphis_a; is++)
        {
            
            ef(is + firsts_a) += weight * (-1.0/dt) * sa * rho_a[0] * phi[0] * phi_ssa(is,0);
            
        }
        
        for (int is = 0; is < nphis_b; is++)
        {
            
            ef(is + firsts_b) += weight * (-1.0/dt) * sb * rho_b[0] * phi[0] * phi_ssb(is,0);
            
        }
        
        return;
    }
    
    
    for (int is = 0; is < nphis_a; is++)
    {
        
        ef(is + firsts_a) += weight * (1.0/dt) * sa * rho_a[0] * phi[0] * phi_ssa(is,0);
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(is + firsts_a, js + firsts_a) += weight * (1.0/dt) * rho_a[0] * phi[0] * phi_ssa(js,0) * phi_ssa(is,0);
        }
        
    }
    
    for (int is = 0; is < nphis_b; is++)
    {
        
        ef(is + firsts_b) += weight * (1.0/dt) * sb * rho_b[0] * phi[0] * phi_ssb(is,0);
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(is + firsts_b, js + firsts_b) += weight * (1.0/dt) * rho_b[0] * phi[0] * phi_ssb(js,0) * phi_ssb(is,0);
        }
        
    }
    
}

void TRMPhaseTransport::Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    
    int nvars = 4; // {p,sa,sb,t}
    
    int sb_a    = 0;
    int sb_b    = 1;
    
    TPZFNMatrix<100,STATE> phi_ssa       = datavec[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb       = datavec[sb_b].phi;
    REAL sa = datavec[sb_a].sol[0][0];
    REAL sb = datavec[sb_b].sol[0][0];
    
    int nphis_a     = phi_ssa.Rows();
    int nphis_b     = phi_ssb.Rows();
    int firsts_a    = 0;
    int firsts_b    = firsts_a + nphis_a;
    
    // Time
    STATE dt = fSimulationData->dt();
    
    // Get the pressure at the integrations points
    long global_point_index = datavec[sb_a].intGlobPtIndex;
    TRMPhaseMemory &point_memory = GetMemory()[global_point_index];
    REAL p_avg_n    = point_memory.p_avg_n();
    REAL p_avg      = point_memory.p_avg();
    
//    REAL sa_avg_n    = point_memory.sa_n();
//    REAL sa_avg      = point_memory.sa();
//    
//    REAL sb_avg_n    = point_memory.sb_n();
//    REAL sb_avg      = point_memory.sb();
    
    //  Average values
    
    //  Computing closure relationship at given average values
    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_avg_n;
    v[1] = sa;
    v[2] = sb;

    
    // Fluid parameters
    TPZManVector<STATE, 10> rho_a,rho_b,l;
    fSimulationData->AlphaProp()->Density(rho_a, v);
    fSimulationData->BetaProp()->Density(rho_b, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi(nvars,0.0);
    //    fSimulationData->Map()->Kappa(datavec[ub].x, K, Kinv, v);
    //    fSimulationData->Map()->phi(datavec[ub].x, phi, v);
    
    // Rock parameters form point memory
    REAL phi_0;
    K       = point_memory.K_0();
    Kinv    = point_memory.Kinv_0();
    phi_0   = point_memory.phi_0();
    phi[0] = phi_0;
    
    
    // Integration point contribution
    
    if(! fSimulationData->IsCurrentStateQ()){
        v[0] = p_avg;
        v[1] = sa;
        v[2] = sb;
        fSimulationData->AlphaProp()->Density(rho_a, v);
        fSimulationData->BetaProp()->Density(rho_b, v);
        fSimulationData->Map()->phi(datavec[sb_a].x, phi, v);
        
        
        for (int is = 0; is < nphis_a; is++)
        {
            
            ef(is + firsts_a) += weight * (-1.0/dt) * sa * rho_a[0] * phi[0] * phi_ssa(is,0);
            
        }
        
        for (int is = 0; is < nphis_b; is++)
        {
            
            ef(is + firsts_b) += weight * (-1.0/dt) * sb * rho_b[0] * phi[0] * phi_ssb(is,0);
            
        }
        
        return;
    }
    
    
    for (int is = 0; is < nphis_a; is++)
    {
        
        ef(is + firsts_a) += weight * (1.0/dt) * sa * rho_a[0] * phi[0] * phi_ssa(is,0);
        
    }
    
    for (int is = 0; is < nphis_b; is++)
    {
        
        ef(is + firsts_b) += weight * (1.0/dt) * sb * rho_b[0] * phi[0] * phi_ssb(is,0);
        
    }
    
}



int TRMPhaseTransport::ClassId() const {
    return -637806378;
}

// -------------------------------------------------------------------------------------------

void TRMPhaseTransport::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMPhaseTransport::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}

// Update element memory by copying the n+1 data to the n data
void TRMPhaseTransport::UpdateMemory()
{
    DebugStop();
//    long nel = fMemory.NElements();
//    for (long el=0; el<nel; el++) {
//        fMemory[el].UpdateSolutionMemory();
//    }
}

