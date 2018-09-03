/**
 * @file
 * @brief TODO:DESCRIBE IT
 * @details BLABLABLA
 *
 * @author Francisco Orlandini
 * @since 2018
 */

#include <TPZSSpStructMatrix.h>
#include <pzsbstrmatrix.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <TPZGmshReader.h>
#include <pzl2projection.h>
#include <pzgeoelrefless.h>
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzelchdiv.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzstepsolver.h"
#include "TPZMatAcousticsTransient.h"
#include "pzgeopoint.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"

void
CreateGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const bool &print,
            const REAL &elSize, const REAL &length, const REAL &height, const std::string &prefix);

void
FilterBoundaryEquations(TPZCompMesh *cmeshHCurl,
                             TPZVec<int64_t> &activeEquations, int &neq,
                             int &neqOriginal);

void
CreateCMesh(TPZCompMesh *&cmesh, TPZGeoMesh *gmesh, int pOrder, const std::string &prefix, const bool &print,
            const TPZVec<int> &matIdVec, const REAL &rho, const REAL &velocity, const REAL &deltaT);

void
RunSimulation(const int &pOrder, const std::string &prefix, const REAL &wZero);


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    const std::string prefix = "results/";//PARAMS
    int pOrder = 2; //PARAMS
    const int nDivIni = 4; //PARAMS
    const int nPcycles = 1;
    const int nHcycles = 1;
    const REAL wZero = 100 * 2 *M_PI;
    boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();
    for (int iP = 0; iP < nPcycles; ++iP, ++pOrder) {
        int nDiv = nDivIni;
        for (int iH = 0; iH < nHcycles; ++iH) {
            std::cout << iH+1<<"/"<<nHcycles<<": Beginning simulation with nEl = "
                      << nDiv * nDiv * 2
                      <<" and p = "<<pOrder<< std::endl;
            RunSimulation(pOrder, prefix,wZero);
            nDiv *= 2;
            std::cout << "************************************" << std::endl;
        }
    }
    boost::posix_time::ptime t2 =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Total time: "<<t2-t1<<std::endl;
    return 0;
}

void RunSimulation(const int &pOrder, const std::string &prefix, const REAL &wZero) {
    // PARAMETROS FISICOS DO PROBLEMA
    const int nThreads = 8; //PARAMS
    const bool l2error = true; //PARAMS
    const bool genVTK = true; //PARAMS
    const bool printG = true;//PARAMS
    const bool printC = true;//PARAMS
    const int postprocessRes = 0;//PARAMS
    REAL rho = 1.3;
    REAL velocity = 340;
    REAL peakTime = 1./100;
    REAL amplitude = 10;
    REAL totalTime = 5 * peakTime;
    REAL elSize = 2 *M_PI*velocity / (10 *wZero),length = 60,height = 8;

    ////////////////////////////////////////////////////////////////////////
    const int64_t nTimeSteps = std::ceil(totalTime/(0.2 *elSize / velocity));
    const REAL deltaT = totalTime/nTimeSteps;
    const REAL cfl = velocity * deltaT/(elSize);
    std::cout<<"CFL: "<<cfl<<std::endl;

    std::cout<<"Creating gmesh... ";
    boost::posix_time::ptime t1_g =
            boost::posix_time::microsec_clock::local_time();
    TPZGeoMesh *gmesh = nullptr;
    std::string meshFileName("wellMesh.geo");
    TPZVec<int> matIdVec;
    CreateGMesh(gmesh, meshFileName, matIdVec, printG, elSize, length,
                height, prefix);

    int matIdSource = 0;
    for (int iMat = 0; iMat< matIdVec.size(); iMat++){
        matIdSource += matIdVec[iMat];
    }
    TPZVec<int> matIdVecCp(matIdVec);
    matIdVec.resize(matIdVec.size()+1);
    matIdVec[0] = matIdSource;
    for (int iMat = 0; iMat< matIdVecCp.size(); iMat++){
        matIdVec[iMat+1] = matIdVecCp[iMat];
    }

    int64_t sourceNodeIndex = -1;
    {
        const REAL sourcePosX = length/2.;
        const REAL sourcePosY = height/2.;
        REAL minDist = 1e6;
        for(int iNode = 0; iNode < gmesh->NNodes(); iNode++){
            const TPZGeoNode & currentNode = gmesh->NodeVec()[iNode];
            const REAL xNode = currentNode.Coord(0);
            const REAL yNode = currentNode.Coord(1);
            const REAL currentDist =
                    (sourcePosX - xNode)*(sourcePosX - xNode) + (sourcePosY - yNode)*(sourcePosY - yNode);
            if(currentDist < minDist){
                minDist = currentDist;
                sourceNodeIndex = currentNode.Id();
            }
        }
        TPZVec<int64_t> nodeIdVec(1,sourceNodeIndex);
        TPZGeoElRefLess<pzgeom::TPZGeoPoint > *zeroDEl =
                new TPZGeoElRefLess<pzgeom::TPZGeoPoint >(nodeIdVec, matIdSource,
                        *gmesh);
        gmesh->BuildConnectivity();
    }

    boost::posix_time::ptime t2_g =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_g-t1_g<<std::endl;

    std::cout<<"Creating cmesh... ";
    boost::posix_time::ptime t1_c =
            boost::posix_time::microsec_clock::local_time();
    TPZCompMesh *cmesh = NULL;
    CreateCMesh(cmesh, gmesh, pOrder, prefix, printC, matIdVec, rho,
                velocity, deltaT);
    boost::posix_time::ptime t2_c =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_c-t1_c<<std::endl;

    TPZAnalysis an(cmesh);
    // configuracoes do objeto de analise
    TPZManVector<int64_t, 1000> activeEquations;
    int neq = 0;
    int neqOriginal = 0;


//#define USING_SKYLINE
#ifdef USING_SKYLINE
    TPZSkylineStructMatrix structMatrix(cmesh);
#else
    TPZSymetricSpStructMatrix structMatrix(cmesh);
#endif
    structMatrix.SetNumThreads(nThreads);
    bool filter = true;
    if(filter){
        FilterBoundaryEquations(cmesh, activeEquations, neq, neqOriginal);
        structMatrix.EquationFilter().SetActiveEquations(activeEquations);
    }else{
        neq = cmesh->NEquations();
    }
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    an.SetStructuralMatrix(structMatrix);

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

    //t=t0
    an.Assemble();

    TPZFMatrix<STATE> currentSol(neq,2,0);
    TPZFMatrix<STATE> scatteredSol(neqOriginal,2,0);
    TPZFMatrix<STATE> allSolutions(neq,nTimeSteps,0);
    auto source = [wZero,peakTime,amplitude](const REAL &time, STATE & val){
        val = (wZero/2*M_PI)*M_SQRT2*
                amplitude*std::exp(0.5 - (wZero/(2*M_PI))*(wZero/(2*M_PI))*M_PI*M_PI*(time - peakTime)*(time - peakTime))*M_PI*(time - peakTime);
    };

    {
        TPZMatAcousticsTransient *mat = dynamic_cast<TPZMatAcousticsTransient *>(cmesh->MaterialVec()[matIdSource]);
        mat->SetSource(source);
    }

    for(int it = 0; it < nTimeSteps; it++){
        for (int imat = 0; imat < matIdVec.NElements(); ++imat) {
            TPZMatAcousticsTransient *mat = dynamic_cast<TPZMatAcousticsTransient *>(cmesh->MaterialVec()[imat]);
            mat->SetCurrentTime(it * deltaT);
        }
        for(int i = 0; i < neq; i++){
            const int prevSol = i - 1 < 0 ? nTimeSteps - i : i-1;
            const int prevPrevSol = i - 2 < 0 ? nTimeSteps - i - 1 : i-2;
            currentSol(i,0) = allSolutions(i,prevPrevSol);
            currentSol(i,1) = allSolutions(i,prevSol);
        }
        if(filter) {
            structMatrix.EquationFilter().Scatter(currentSol, scatteredSol);
            an.LoadSolution(scatteredSol);
        }else{
            an.LoadSolution(currentSol);
        }
        an.AssembleResidual();
        an.Solve();
        TPZFMatrix<STATE>& stepSol = an.Solution();
        for(int i = 0; i < neq; i++){
            allSolutions(i,it) = stepSol(i,0);
        }
    }

    if (genVTK) {
        std::cout<<"Post processing... "<<std::endl;

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        std::string plotfile = prefix+"sol";
        plotfile.append(std::to_string(cmesh->NElements()));
        plotfile.append(".vtk");

        an.DefineGraphMesh(2, scalnames, vecnames,
                           plotfile);  // define malha grafica
        int postProcessResolution = postprocessRes; // define resolucao do pos processamento

        TPZFMatrix<STATE> currentSol(neq,1);
        TPZFMatrix<STATE> scatteredSol(neqOriginal,1);
        for(int iTime = 0; iTime < nTimeSteps; iTime++){
            std::cout<<"\rtime: "<<iTime+1<<" out of "<<nTimeSteps<<std::flush;
            for(int iPt = 0; iPt < neq; iPt++){
                currentSol(iPt,0) = allSolutions(iPt,iTime);
            }
            if(filter){
                structMatrix.EquationFilter().Scatter(currentSol, scatteredSol);
                an.LoadSolution(scatteredSol);
            }else{
                an.LoadSolution(currentSol);
            }
            an.PostProcess(postProcessResolution);
        }
        std::cout<<std::endl<<" Done!"<<std::endl;
    }
    gmesh->SetReference(nullptr);
    for(int icon = 0; icon < cmesh->NConnects(); icon++){
        TPZConnect &con = cmesh->ConnectVec()[icon];
        con.RemoveDepend();
    }
    cmesh->SetReference(nullptr);
    delete cmesh;
    cmesh=nullptr;
    delete gmesh;
    gmesh=nullptr;
}

void
CreateGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const bool &print,
            const REAL &elSize, const REAL &length, const REAL &height, const std::string &prefix) {

    std::ostringstream str_elSize;
    str_elSize << std::setprecision(20) << elSize;
    std::ostringstream str_length;
    str_length<< std::setprecision(20) << length;
    std::ostringstream str_height;
    str_height<< std::setprecision(20) << height;

    std::string command = "gmsh " + mshFileName + " -2 -match ";
    command += " -v 3 ";
    command += " -setnumber length "+str_length.str();
    command += " -setnumber height "+str_height.str();
    command += " -setnumber el_size "+str_elSize.str();
    command += " -o " + prefix + "wellMesh.msh";
    std::cout<<"Generating mesh with: "<<std::endl<<command<<std::endl;

    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), 128, pipe.get()) != nullptr){
        result += buffer.data();
    }
    std::cout<<result<<std::endl;


    TPZGmshReader meshReader;
    gmesh = meshReader.GeometricGmshMesh(prefix+"wellMesh.msh");
#ifdef PZDEBUG
	TPZCheckGeom * Geometrytest = new TPZCheckGeom(gmesh);
	int isBadMeshQ = Geometrytest->PerformCheck();

	if (isBadMeshQ) {
		DebugStop();
	}
#endif
    auto matIds = meshReader.fMaterialDataVec;
    matIdVec.Resize(2);
    //fPZMaterialId[dimension][name]
    matIdVec[0] = meshReader.fPZMaterialId[2]["water"];
    matIdVec[1] = meshReader.fPZMaterialId[1]["bound"];

    if(print){
        std::string meshFileName = prefix + "gmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        gmesh->Print(outTXT);
        outTXT.close();
        outVTK.close();
    }

    return;
}

void FilterBoundaryEquations(TPZCompMesh *cmeshHCurl,
                             TPZVec<int64_t> &activeEquations, int &neq,
                             int &neqOriginal) {
    TPZManVector<int64_t, 1000> allConnects;
    std::set<int64_t> boundConnects;

    for (int iel = 0; iel < cmeshHCurl->NElements(); iel++) {
        TPZCompEl *cel = cmeshHCurl->ElementVec()[iel];
        if (cel == NULL) {
            continue;
        }
        if (cel->Reference() == NULL) {
            continue;
        }
        TPZBndCond *mat = dynamic_cast<TPZBndCond *>(cmeshHCurl->MaterialVec()[cel->Reference()->MaterialId()]);
        if (mat && mat->Type() == 0) {
            std::set<int64_t> boundConnectsEl;
            cel->BuildConnectList(boundConnectsEl);

            for (std::set<int64_t>::iterator iT = boundConnectsEl.begin();
                 iT != boundConnectsEl.end(); iT++) {
                const int64_t val = *iT;
                if (boundConnects.find(val) == boundConnects.end()) {
                    boundConnects.insert(val);
                }
            }
        }
    }
    neqOriginal = cmeshHCurl->NEquations();
    for (int iCon = 0; iCon < cmeshHCurl->NConnects(); iCon++) {
        if (boundConnects.find(iCon) == boundConnects.end()) {
            TPZConnect &con = cmeshHCurl->ConnectVec()[iCon];
            int seqnum = con.SequenceNumber();
            int pos = cmeshHCurl->Block().Position(seqnum);
            int blocksize = cmeshHCurl->Block().Size(seqnum);
            if (blocksize == 0)
                continue;

            int vs = activeEquations.size();
            activeEquations.Resize(vs + blocksize);
            for (int ieq = 0; ieq < blocksize; ieq++) {
                activeEquations[vs + ieq] = pos + ieq;
                neq++;
            }
        }
    }
    neqOriginal = cmeshHCurl->NEquations();
    std::cout << "# equations(before): " << neqOriginal << std::endl;
    std::cout << "# equations(after): " << neq << std::endl;
    return;
}

void
CreateCMesh(TPZCompMesh *&cmesh, TPZGeoMesh *gmesh, int pOrder, const std::string &prefix, const bool &print,
            const TPZVec<int> &matIdVec, const REAL &rho, const REAL &velocity, const REAL &deltaT) {
    const int dim = 2;   // dimensao do problema
    const int matIdSource = matIdVec[0]; // define id para um material(formulacao fraca)
    const int matId = matIdVec[1]; // define id para um material(formulacao fraca)
    const int bc0 = matIdVec[matIdVec.size()-1];  // define id para um material(cond contorno dirichlet)
    enum {
        dirichlet = 0,
        neumann,
        mixed
    }; // tipo da condicao de contorno do problema
    // Criando material

    cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha
    TPZMatAcousticsTransient *matAcoustics = NULL;
    //source
    auto fakeSource = [](const REAL&time, STATE&val){
        val = 0;
    };
    matAcoustics = new TPZMatAcousticsTransient(matIdVec[0], rho, velocity);
    matAcoustics->SetDeltaT(deltaT);
    matAcoustics->SetSource(fakeSource);
    cmesh->InsertMaterialObject(matAcoustics);
    //water
    matAcoustics = new TPZMatAcousticsTransient(matIdVec[1], rho, velocity);
    matAcoustics->SetDeltaT(deltaT);
    matAcoustics->SetSource(fakeSource);
    cmesh->InsertMaterialObject(matAcoustics);

    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *bCondDir = matAcoustics->CreateBC(
            matAcoustics, bc0, dirichlet, val1, val2); // cria material que implementa a
    // condicao de contorno de
    // dirichlet

    cmesh->InsertMaterialObject(bCondDir); // insere material na malha
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    if(print){
        std::string meshFileName = prefix + "cmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh,outVTK,true);
        cmesh->Print(outTXT);
        outTXT.close();
        outVTK.close();
    }
}
