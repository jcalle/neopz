#include <map>
#include "pzvec.h"
#include "TPZMatAcousticsFourier.h"
#include "TPZAcousticGeoMesher.h"
#include "pzcmesh.h"

class TPZGeoMesh;
/**
 * A class that creates a computational mesh based on a previously generated geometric mesh.
 * The geometric mesh must have been generated using a TPZAcousticGeoMesher in order to assure that all
 * relevant data structures have been properly filled.
 *
 */
class TPZAcousticCompMesher{
public:

    /**
     * This constructor will not be generated
     */
    TPZAcousticCompMesher() = delete;

    TPZAcousticCompMesher(TPZAcousticGeoMesher * geoMesh, bool isAxisymmetric);

    void CreateFourierMesh(const int &porder);

    void CreateTransientMesh(const int &porder);

    /**
     * Prints the computational mesh in both .vtk and .txt formats.
     * @param gmesh mesh to be printed (in)
     * @param fileName where to print the mesh (extensions will be added to fileName)
     * @param prefix directory to output the printed mesh (in)
     * @param printVTK whether to print the mesh in .vtk format (in)
     * @param printTxt whether to print the mesh in .txt format (in)
     */
    void PrintMesh(const std::string &fileName, const std::string &prefix, bool printVTK, bool printTxt);
protected:
    template<class T>
    void CreateCompMesh(const int &porder);

    void FilterBoundaryEquations(TPZVec<int64_t> &activeEquations, int64_t &neq, int64_t &neqOriginal) ;

    TPZCompMesh *fCmesh;

    const TPZAcousticGeoMesher *fGeoMesh;

    const bool fIsAxisymetric;

    TPZVec<int64_t > fActiveEquations;

    int64_t fNeqReduced;

    int64_t fNeqOriginal;
};

template<class T>
void TPZAcousticCompMesher::CreateCompMesh(const int & pOrder) {
#ifdef PZDEBUG
    if(fCmesh!=nullptr){
        PZError<<"You are trying to create a computational mesh over an existent one. Aborting now..."<<std::endl;
        DebugStop();
    }
#endif
    TPZManVector<int, 8> sourceMatIdVec(1);//TODO:this is wrong
    std::map<int,REAL> rhoMap = fGeoMesh->GetDensityMap();
    std::map<int,REAL> velocityMap = fGeoMesh->GetVelocityMap();
    TPZManVector<int, 8> matIdVec(rhoMap.size());
    int i = 0;
    for (auto itMap : rhoMap){
        matIdVec[i] = itMap.first;
        i++;
    }
    TPZManVector<int, 8> boundTypeVec(1,0);//TODO:this is wrong
    for (int i = 0; i < 1; i++) {
        sourceMatIdVec[i] = matIdVec[i];
    }

    TPZManVector<int, 8> volMatIdVec(matIdVec.size() - boundTypeVec.size() - sourceMatIdVec.size());
    for (int i = 0; i < volMatIdVec.size(); i++) {
        volMatIdVec[i] = matIdVec[i + sourceMatIdVec.size()];
    }
    TPZManVector<int, 8> boundMatIdVec(boundTypeVec.size());
    for (int i = 0; i < boundMatIdVec.size(); i++) {
        boundMatIdVec[i] = matIdVec[volMatIdVec.size() + sourceMatIdVec.size()  + i];
    }


//    const int64_t outerMaterialIndex = volMatIdVec[volMatIdVec.size() - 1];
    //TODO: Introduce axisymmetry
    const int dim = 2;   // dimensao do problema
    TPZCompMesh * cmesh = new TPZCompMesh(fGeoMesh->GetMesh());
    cmesh->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha
    T *matAcoustics = nullptr;//, *matOuter = nullptr;

    for (int i = 0; i < volMatIdVec.size(); ++i) {
        const REAL rho = rhoMap.at(volMatIdVec[i]);
        const REAL velocity = velocityMap.at(volMatIdVec[i]);
        matAcoustics = new T(volMatIdVec[i], rho, velocity);
        cmesh->InsertMaterialObject(matAcoustics);
//        if(volMatIdVec[i] == outerMaterialIndex){
//            rhoOuter = rho;
//            velocityOuter = velocity;
//            matOuter = matAcoustics;
//        }
    }

    for (int i = 0; i < sourceMatIdVec.size(); ++i) {
        const REAL rho = rhoMap.at(volMatIdVec[i]);
        const REAL velocity = velocityMap.at(volMatIdVec[i]);
        matAcoustics = new T(sourceMatIdVec[i], rho, velocity);//rho and v wont be used
        cmesh->InsertMaterialObject(matAcoustics);
    }

    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *bcond = nullptr;
    for (int i = 0; i < boundMatIdVec.size(); i++) {
        bcond = matAcoustics->CreateBC(matAcoustics, boundMatIdVec[i], boundTypeVec[i], val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    FilterBoundaryEquations(fActiveEquations, fNeqReduced, fNeqOriginal);
}
