#include "TPZAcousticTimeDomainAnalysis.h"
#include "TPZAcousticCompMesher.h"
#include "TPZSpStructMatrix.h"
#include "TPZMatAcousticsTransient.h"
#include <pzgeoelrefless.h>
#include <pzgeopoint.h>
#include "pzbndcond.h"

#include <boost/date_time/posix_time/posix_time.hpp>

TPZAcousticTimeDomainAnalysis::TPZAcousticTimeDomainAnalysis(TPZAcousticCompMesher *compMesher, const int &nThreads,
const REAL &deltaT, const int &nTimeSteps, const bool &filter)
        : TPZAcousticAnalysis(compMesher, nThreads, deltaT, nTimeSteps, filter) {
#ifdef PZDEBUG
    if(fCompMesher->fMeshType != TPZAcousticCompMesher::ECompMeshTypes::timeDomain){
        PZError<<"It seems you are trying to run a time domain problem with a frequencydomain computational mesh.";
        PZError<<"Aborting..."<<std::endl;
    }
#endif
    auto dummySource = [](const REAL &time, STATE &val){
        val = 0;
    };
    for (auto itMat : compMesher->fCmesh->MaterialVec()){
        TPZMatAcousticsTransient *mat = dynamic_cast<TPZMatAcousticsTransient *>(itMat.second);
        if(mat){
            mat->SetSource(dummySource);
        }
    }
    fAmIready = false;
}


void TPZAcousticTimeDomainAnalysis::InitializeComputations() {

    this->InitializeSolver();

    TPZCompMesh *cmesh = fCompMesher->fCmesh;

    for(auto itMat : cmesh->MaterialVec()){
        TPZMatAcousticsTransient *mat = dynamic_cast<TPZMatAcousticsTransient *>(itMat.second);
        if(mat){
            mat->SetDeltaT(fDeltaT);
        }
    }

    fPzAnalysis.Assemble();

    fAmIready = true;
}

void TPZAcousticTimeDomainAnalysis::RunSimulationSteps() {

    if(! fAmIready ){
        PZError<<"TPZAcousticTimeDomainAnalysis::InitializeComputations should be called before RunSimulationSteps."<<std::endl;
        DebugStop();
    }
    TPZCompMesh *cmesh = fCompMesher->fCmesh;
    TPZGeoMesh *gmesh = fCompMesher->fGeoMesh->fGmesh;

    const int matIdSource = fCompMesher->fGeoMesh->fMatIdSource;

    TPZMatAcousticsTransient *mat = dynamic_cast<TPZMatAcousticsTransient *>(cmesh->MaterialVec()[matIdSource]);
    mat->SetSource(fTimeDomainSource);

    TPZCompEl * sourceCompel = nullptr;
    {
        TPZGeoElRefLess<pzgeom::TPZGeoPoint> *sourceEl = nullptr;
        for (int iEl = 0; iEl < gmesh->NElements(); ++iEl) {
            sourceEl = dynamic_cast<TPZGeoElRefLess<pzgeom::TPZGeoPoint>*>(gmesh->Element(iEl));
            if(sourceEl){
                break;
            }
        }
        if(sourceEl == nullptr){
            std::cout<<"could not find source element"<<std::endl;
            DebugStop();
        }
        sourceCompel = sourceEl->Reference();
    }


    uint solSize = fNeqOriginal;
    TPZFMatrix<STATE> currentSol(solSize,2,0);
    fProbeSolution.Resize(fNTimeSteps,2);
    fTimeDomainSolution.Resize(solSize,fNTimeSteps);
    for(int it = 0; it < fNTimeSteps; it++){
        std::cout<<"time step "<<it+1<<" out of "<<fNTimeSteps <<"\t";
        for(auto imat : cmesh->MaterialVec()) {
            auto *mat = dynamic_cast<TPZMatAcousticsTransient *>(imat.second);
            auto *bound = dynamic_cast<TPZBndCond *>(imat.second);
            if(mat!=nullptr){
                mat->SetCurrentTime(it*fDeltaT);
            }
            else if(bound== nullptr){
                DebugStop();
            }
        }
        const int prevSol = it - 1 < 0 ? fNTimeSteps - it - 1: it-1;
        const int prevPrevSol = it - 2 < 0 ? fNTimeSteps - it - 2 : it-2;

        for(int i = 0; i < solSize; i++){
            currentSol(i,0) = fTimeDomainSolution(i,prevSol);
            currentSol(i,1) = fTimeDomainSolution(i,prevPrevSol);
        }
        fPzAnalysis.LoadSolution(currentSol);
        fPzAnalysis.AssembleResidual();
        fPzAnalysis.Solve();
        TPZFMatrix<STATE>& stepSol = fPzAnalysis.Solution();
        std::cout<<std::endl;
        for(int i = 0; i < solSize; i++){
//            std::cout<<"stepSol("<<i<<",0)"<<stepSol(i,0)<<std::endl;
            fTimeDomainSolution(i,it) = stepSol(i,0);
        }

        //get probe solution
        fProbeSolution(it,0) = it*fDeltaT;
        TPZVec<REAL> qsi(1,0);
        int var = 1;
        TPZVec<STATE> sol(1,0);
        sourceCompel->Solution(qsi,var,sol);
        fProbeSolution(it,1) = sol[0];
    }
}

void TPZAcousticTimeDomainAnalysis::SetUpGaussianSource(const REAL &wZero, const REAL &peakTime, const REAL &amplitude) {
#ifdef PZDEBUG
    if(fCompMesher->fMeshType != TPZAcousticCompMesher::ECompMeshTypes::timeDomain){
        PZError<<"Trying to define a frequency domain source for a problem that is not on time domain."<<std::endl;
        DebugStop();
    }
    const REAL tol = 1e-12;
    if(wZero < tol || peakTime < tol || amplitude < tol){
        PZError<<"The time domain gaussian source is not properly defined"<<std::endl;
        DebugStop();
    }
#endif
    auto source = [wZero,peakTime,amplitude](const REAL &time, STATE & val){
        val = 0.5 * amplitude * (peakTime-time)*wZero*wZero*std::exp(-.25 * (peakTime-time)*(peakTime-time)*wZero*wZero);
//        val = (wZero/2*M_PI)*M_SQRT2*                amplitude*std::exp(0.5 - (wZero/(2*M_PI))*(wZero/(2*M_PI))*M_PI*M_PI*(time - peakTime)*(time - peakTime))*M_PI*(time - peakTime);
    };

    fTimeDomainSource = source;
}


void TPZAcousticTimeDomainAnalysis::InitializeSolver() {

#ifdef STATE_COMPLEX
    fStepSolver.SetDirect(ELU);
#else
    fStepSolver.SetDirect(ELDLt);
#endif
    fPzAnalysis.SetSolver(fStepSolver);
}