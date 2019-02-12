#include "TPZAcousticFreqDomainAnalysis.h"
#include "TPZAcousticCompMesher.h"
#include "TPZSpStructMatrix.h"
#include <pzgeoelrefless.h>
#include <pzgeopoint.h>
#include "pzbndcond.h"

#include <boost/date_time/posix_time/posix_time.hpp>

TPZAcousticFreqDomainAnalysis::TPZAcousticFreqDomainAnalysis(TPZAcousticCompMesher *compMesher, const int &nThreads,
const REAL &deltaT, const int &nTimeSteps, const bool &filter)
        : TPZAcousticAnalysis(compMesher, nThreads, deltaT, nTimeSteps, filter) {
#ifdef PZDEBUG
    if(fCompMesher->fMeshType != TPZAcousticCompMesher::ECompMeshTypes::freqDomain){
        PZError<<"It seems you are trying to run a frequency domain problem with a time domain computational mesh.";
        PZError<<"Aborting..."<<std::endl;
    }
#endif
    fAFactor = -1;
    fWMax = -1;
    fNSamples = -1;
    fAmIready = false;
}


void TPZAcousticFreqDomainAnalysis::InitializeComputations() {

    this->InitializeSolver();

    TPZCompMesh *cmesh = fCompMesher->fCmesh;

    std::set<int> matsNonPml;
    std::set<int> matSource;

    int matIdSource = fCompMesher->fGeoMesh->fMatIdSource;
    for (auto itMap : cmesh->MaterialVec()) {
        TPZMatAcousticsFourier *matH1 = dynamic_cast<TPZMatAcousticsFourier *>(itMap.second);
        if (matH1 != nullptr) {
            if (matH1->Id() == matIdSource) {
                matSource.insert(matH1->Id());
                continue;
            }
            matsNonPml.insert(matH1->Id());
        } else {
            TPZBndCond *matBound = dynamic_cast<TPZBndCond *>(itMap.second);
            if( matBound == nullptr){
                PZError << "Something is wrong... " << __PRETTY_FUNCTION__ << " could not identify material with id ";
                PZError << itMap.first << std::endl;
            }
        }
    }


    fMatrices.Resize(2);//M, K
    //GET M MATRIX
    {
        TPZAutoPointer<TPZStructMatrix> structMatrixPtr = fPzAnalysis.StructMatrix();
        structMatrixPtr->SetMaterialIds(matsNonPml);
        //set structMatrix M
        for (std::set<int>::iterator it = matsNonPml.begin(); it != matsNonPml.end(); ++it) {
            TPZMatAcousticsFourier *matH1 = dynamic_cast<TPZMatAcousticsFourier *>(cmesh->MaterialVec()[*it]);
            matH1->SetAssemblingMatrix(TPZMatAcousticsFourier::M);
        }
        fPzAnalysis.Assemble();
        //get structMatrix M
        auto mat = dynamic_cast<TPZFYsmpMatrix<STATE> *>( fPzAnalysis.Solver().Matrix().operator->());
        if (!mat) {
            DebugStop();
        }
        fMatrices[0] = *mat;
    }

    //GET K MATRIX
    {
        TPZAutoPointer<TPZStructMatrix> structMatrixPtr = fPzAnalysis.StructMatrix();
        structMatrixPtr->SetMaterialIds(matsNonPml);
        //set structMatrix M
        for (std::set<int>::iterator it = matsNonPml.begin(); it != matsNonPml.end(); ++it) {
            TPZMatAcousticsFourier *matH1 = dynamic_cast<TPZMatAcousticsFourier *>(cmesh->MaterialVec()[*it]);
            matH1->SetAssemblingMatrix(TPZMatAcousticsFourier::K);
        }
        fPzAnalysis.Assemble();
        //get structMatrix M
        auto mat = dynamic_cast<TPZFYsmpMatrix<STATE> *>( fPzAnalysis.Solver().Matrix().operator->());
        if (!mat) {
            DebugStop();
        }
        fMatrices[1] = *mat;
    }
}

void TPZAcousticFreqDomainAnalysis::SetUpFourierSettings(const REAL &wMax, const int &nSamples,
                                                         const REAL &aFactorUser) {
    if(nSamples < 2){
        std::cout<<"The number of frequency samples does not seem right, aborting..."<<std::endl;
        DebugStop();
    }


    const REAL wSample = wMax/nSamples;
    fWMax = wMax;
    fNSamples = nSamples;
    if(aFactorUser < 0){
        std::cout<<"aFactor = "<<aFactorUser<<" using default value..."<<std::endl;
        fAFactor = log(50)/(2*M_PI/wSample);
    }else{
        fAFactor = aFactorUser;
    }
    fAmIready = true;
}
void TPZAcousticFreqDomainAnalysis::RunSimulationSteps() {

    if(! fAmIready ){
        PZError<<"TPZAcousticFreqDomainAnalysis::SetUpFourierSettings should be called before RunSimulationSteps."<<std::endl;
        DebugStop();
    }
    int &nSamples = fNSamples;
    REAL &wMax = fWMax;
    const REAL wSample = wMax/nSamples;
    TPZVec<int64_t> ia;
    TPZVec<int64_t> ja;
    TPZVec<STATE> aM, aK, aFinal;
    TPZFYsmpMatrix<STATE> &matM = fMatrices[0];
    TPZFYsmpMatrix<STATE> &matK = fMatrices[1];
    matM.GetData(ia,ja,aM);
    matK.GetData(ia,ja,aK);
    aFinal.Resize(aK.size());



    ///////
    TPZCompEl * sourceCompel = nullptr;
    TPZGeoMesh *gmesh = fCompMesher->fGeoMesh->fGmesh;
    TPZCompMesh *cmesh = fCompMesher->fCmesh;
    const int &matIdSource = fCompMesher->fGeoMesh->fMatIdSource;
    std::set<int>matSource;
    matSource.insert(matIdSource);
    TPZMatAcousticsFourier *matSourcePtr = dynamic_cast<TPZMatAcousticsFourier *>(cmesh->MaterialVec()[matIdSource]);
    if(matSourcePtr == nullptr){
        DebugStop();
    }
    matSourcePtr->SetSource(fFreqDomainSource);
    {
        TPZGeoElRefLess<pzgeom::TPZGeoPoint> *sourceEl = nullptr;
        for (int iEl = 0; iEl < gmesh->NElements(); ++iEl) {
            if (gmesh->Element(iEl)->MaterialId() == matIdSource){
                sourceEl = dynamic_cast<TPZGeoElRefLess<pzgeom::TPZGeoPoint>*>(gmesh->Element(iEl));
                if(sourceEl){
                    break;
                }else{
                    DebugStop();
                }
            }
        }
        if(sourceEl == nullptr){
            std::cout<<"could not find source element"<<std::endl;
            DebugStop();
        }
        sourceCompel = sourceEl->Reference();
    }

    fProbeSolution.Resize(nSamples-1,2);


//    const int solSize = fFilterBoundaryEquations ? fCompMesher->fNeqReduced : fCompMesher->fNeqOriginal;
    const int solSize = fNeqOriginal;
    fTimeDomainSolution.Resize(solSize,fNTimeSteps);
    for(int iW = 0; iW < nSamples-1; iW++){

        TPZAutoPointer<TPZMatrix<STATE>> matFinal(new TPZFYsmpMatrix<STATE>(matK.Rows(),matK.Cols()));

        const STATE currentW = (iW+1) * wSample + SPZAlwaysComplex<STATE>::type(0,1)*fAFactor;
        //-w^2 M + K
        boost::posix_time::ptime t1_sum =
                boost::posix_time::microsec_clock::local_time();
        //Calculate final Matrix
        TPZFYsmpMatrix<STATE> * matFinalDummy = dynamic_cast<TPZFYsmpMatrix<STATE> *>(matFinal.operator->());
        for(int i = 0; i < aFinal.size(); i++){
            aFinal[i] = -1.*currentW*currentW*aM[i] + aK[i];
        }
        matFinalDummy->SetData(ia,ja,aFinal);


        boost::posix_time::ptime t2_sum =
                boost::posix_time::microsec_clock::local_time();
        std::cout<<"Summing matrices took "<<t2_sum-t1_sum<<std::endl;


        //pml elements
        TPZAutoPointer<TPZStructMatrix> structMatrixPtr = fPzAnalysis.StructMatrix();
        //assemble pml
//        if(nPmls > 0){
//            structMatrixPtr->SetMaterialIds(matsPml);
//            for(auto itMap : matsPml) {
//                TPZMatAcousticsPml *mat = dynamic_cast<TPZMatAcousticsPml *>(cmesh->FindMaterial(itMap));
//                if (mat != nullptr) {
//                    mat->SetW(currentW);
//                } else {
//                    DebugStop();
//                }
//            }
//            structMatrixPtr->Assemble(matFinal,rhsFake,nullptr);
//        }


        structMatrixPtr->SetMaterialIds(matSource);
        matSourcePtr->SetW(currentW);//FONTE IGUAL A DO COMSOL
        //assemble load vector

        fPzAnalysis.AssembleResidual();
        //solve system
        fPzAnalysis.Solver().SetMatrix(matFinal);

        std::cout<<"Beginning solver step "<<iW+1<<" out of "<<nSamples - 1<<std::endl;
        boost::posix_time::ptime t1_solve =
                boost::posix_time::microsec_clock::local_time();
        fPzAnalysis.Solve();
        boost::posix_time::ptime t2_solve =
                boost::posix_time::microsec_clock::local_time();
        std::cout<<"solver step "<<iW+1<<" out of "<<nSamples - 1<<" took "<<t2_solve-t1_solve<<std::endl;
        //get solution
        TPZFMatrix<STATE> &currentSol = fPzAnalysis.Solution();//p omega

        //get spectrum
        fProbeSolution(iW,0) = currentW;
        TPZVec<REAL> qsi(1,0);
        int var = 1;
        TPZVec<STATE> sol(1,0);
        sourceCompel->Solution(qsi,var,sol);
        fProbeSolution(iW,1) = sol[0];
        //inverse transform


        for(int iPt = 0; iPt < solSize; iPt++){
            for(int iTime = 0; iTime < fNTimeSteps; iTime++){
                const REAL currentTime = (REAL)iTime * fDeltaT;
                fTimeDomainSolution(iPt,iTime) +=
                        std::exp(fAFactor*currentTime)*//fator para compensar o shift na frequência
                        wSample *//fator para compensar o trem de impulso com intensidade wSample
                        0.5*(1.+std::cos(std::real(currentW)*M_PI/wMax))*//Hanning windowing, atenua frequencias altas
                        std::real(currentSol(iPt,0)) * std::cos(-1. * currentTime * std::real(currentW))*
                        M_2_SQRTPI * M_SQRT2 * M_SQRT1_2;//1/sqrt(2pi), para ser coerente com a
                // transformada feita no mathematica
            }
        }
    }
}

void TPZAcousticFreqDomainAnalysis::SetUpGaussianSource(const REAL &wZero, const REAL &peakTime, const REAL &amplitude) {
#ifdef PZDEBUG
    if(fCompMesher->fMeshType != TPZAcousticCompMesher::ECompMeshTypes::freqDomain){
        PZError<<"Trying to define a frequency domain source for a problem that is not on frequency domain."<<std::endl;
        DebugStop();
    }
    const REAL tol = 1e-12;
    if(wZero < tol || peakTime < tol || amplitude < tol){
        PZError<<"The frequency domain gaussian source is not properly defined"<<std::endl;
        DebugStop();
    }
#endif
    auto source = [wZero,peakTime,amplitude](const STATE &currentW, STATE & val){
        val = 0.;
        val += -1. * SPZAlwaysComplex<STATE>::type(0,1)  * amplitude * currentW * M_SQRT2/wZero;
        val*= exp(SPZAlwaysComplex<STATE>::type(0,1)  * currentW * peakTime - (currentW/wZero)*(currentW/wZero));
    };

    fFreqDomainSource = source;
}


void TPZAcousticFreqDomainAnalysis::InitializeSolver() {
    fStepSolver.SetDirect(ELU);
    fPzAnalysis.SetSolver(fStepSolver);
}