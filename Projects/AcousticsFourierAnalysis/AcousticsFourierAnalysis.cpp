/**
 * @file
 * @brief TODO:DESCRIBE IT
 * @details BLABLABLA
 *
 * @author Francisco Orlandini
 * @since 2018
 */

#include <TPZSSpStructMatrix.h>
#include <TPZSpStructMatrix.h>
#include <pzskylstrmatrix.h>
#include <pzsbstrmatrix.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <TPZGmshReader.h>
//#include "TPZMatHelmholtz2D.h"
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
#include "TPZMatAcousticsH1.h"
#include "TPZMatAcousticsPml.h"

void
CreateGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const bool &print,
            const REAL &elSize, const REAL &length, const REAL &height, const REAL &pmlLength, const std::string &prefix);
void loadVec(const TPZVec<REAL> &coord, TPZVec<STATE> &val) {
    val.Resize(3, 0.);
    //  val[0] = 3 - coord[1] * coord[1];
    //  val[1] = 3 - coord[0] * coord[0];
    val[0] = (2 * M_PI * M_PI + 1.) * M_PI *
             cos(M_PI * coord[0]) * sin(M_PI * coord[1]);
    val[1] = (2 * M_PI * M_PI + 1.) * M_PI * (-1.) *
             sin(M_PI * coord[0]) * cos(M_PI * coord[1]);
}

//void exactSol(const TPZVec<REAL> &coord, TPZVec<STATE> &result,
//              TPZFMatrix<STATE> &curl) {
//    result.Resize(3, 0.);
//    result[0] = M_PI * cos(M_PI * coord[0]) * sin(M_PI * coord[1]);
//    result[1] = M_PI * (-1.) * sin(M_PI * coord[0]) * cos(M_PI * coord[1]);
//
//    curl.Resize(1, 1);
//    curl(0, 0) = -2 * M_PI * M_PI * cos(M_PI * coord[0]) * cos(M_PI * coord[1]);
//}

void FilterBoundaryEquations(TPZCompMesh *cmeshHCurl,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);

void
CreateCMesh(TPZCompMesh *&cmesh, TPZGeoMesh *gmesh, int pOrder, void (&loadVec)(const TPZVec<REAL> &, TPZVec<STATE> &),
            const REAL &alphaPml, const std::string &prefix, const bool &print,
            const TPZVec<int> &matIdVec, const REAL &rho, const REAL &velocity);

void RunSimulation(const int &nDiv, const int &pOrder, const std::string &prefix, const REAL &wZero);


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    const std::string prefix = "results/";//PARAMS
    int pOrder = 1; //PARAMS
    const int nDivIni = 4; //PARAMS
    const int nPcycles = 4;
    const int nHcycles = 7;
    const REAL wZero = 100 * 2 *M_PI;
    boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();
    for (int iP = 0; iP < nPcycles; ++iP, ++pOrder) {
        int nDiv = nDivIni;
        for (int iH = 0; iH < nHcycles; ++iH) {
            std::cout << iH+1<<"/"<<nHcycles<<": Beginning simulation with nEl = "
                      << nDiv * nDiv * 2
                      <<" and p = "<<pOrder<< std::endl;
            RunSimulation(nDiv, pOrder, prefix,wZero);
            nDiv *= 2;
            std::cout << "************************************" << std::endl;
        }
    }
    boost::posix_time::ptime t2 =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Total time: "<<t2-t1<<std::endl;
    return 0;
}

void RunSimulation(const int &nDiv, const int &pOrder, const std::string &prefix, const REAL wZero) {
    // PARAMETROS FISICOS DO PROBLEMA
    const int nThreads = 8; //PARAMS
    const bool l2error = true; //PARAMS
    const bool genVTK = false; //PARAMS
    const bool printG = false;//PARAMS
    const bool printC = false;//PARAMS
    const int postprocessRes = 0;//PARAMS
    REAL elSize = 4,length = 60,height = 20,pmlLength = 20;
    REAL alphaPML;
    REAL rho;
    REAL velocity;

    std::cout<<"Creating gmesh... ";
    boost::posix_time::ptime t1_g =
            boost::posix_time::microsec_clock::local_time();
    TPZGeoMesh *gmesh = NULL;
    std::string meshFileName;
    TPZVec<int> matIdVec;
    CreateGMesh(gmesh, meshFileName, matIdVec, printG, elSize, length,
                height, pmlLength, prefix);
    boost::posix_time::ptime t2_g =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_g-t1_g<<std::endl;

    std::cout<<"Creating cmesh... ";
    boost::posix_time::ptime t1_c =
            boost::posix_time::microsec_clock::local_time();
    TPZCompMesh *cmesh = NULL;
    CreateCMesh(cmesh, gmesh, pOrder, loadVec, alphaPML, prefix, printC, matIdVec, rho, velocity);
    boost::posix_time::ptime t2_c =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_c-t1_c<<std::endl;

    TPZAnalysis an(cmesh);
    // configuracoes do objeto de analise
    TPZManVector<long, 1000> activeEquations;
    int neq = 0;
    int neqOriginal = 0;

//  TPZAutoPointer<TPZStructMatrix> strmtrx;
//  strmtrx = new TPZSpStructMatrix(cmeshHCurl);
//  strmtrx->SetNumThreads(nThreads);
//  FilterBoundaryEquations(cmeshHCurl, activeEquations, neq, neqOriginal);
//  strmtrx->EquationFilter().SetActiveEquations(activeEquations);
//  an.SetStructuralMatrix(strmtrx);

    TPZSymetricSpStructMatrix matrix(cmesh);
    //TPZSBandStructMatrix matrix(cmeshHCurl);
    matrix.SetNumThreads(nThreads);
    FilterBoundaryEquations(cmesh, activeEquations, neq, neqOriginal);
    matrix.EquationFilter().SetActiveEquations(activeEquations);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.SetStructuralMatrix(matrix);
    ////////////////////////////////////////////////////////////////////////
    const REAL wMax = 3*wZero;
    const int nSamples = 100;
    const REAL wSample = wMax/nSamples;

    ////////////////////////////////////////////////////////////////////////
    std::set<int> matsPml;
    std::set<int> matsNonPml;
    for(auto itMap : cmesh->MaterialVec()){
        TPZMatAcousticsPml *matPml = dynamic_cast<TPZMatAcousticsPml *>(itMap.second);
        TPZMatAcousticsH1 *matH1 = dynamic_cast<TPZMatAcousticsH1 *>(itMap.second);
        if(matH1!=nullptr && matPml==nullptr){
            matsNonPml.insert(matH1->Id());
        }else if(matPml !=nullptr){
            matsPml.insert(matPml->Id());
        }
    }

    {
        TPZAutoPointer<TPZStructMatrix> structMatrixPtr = an.StructMatrix();
        structMatrixPtr->SetMaterialIds(matsNonPml);
        //set matrix M
        an.Assemble();
        //get matrix M

        //set matrix K
        an.Assemble();
        //get matrix K
    }



    ////////////////////////////////////////////////////////////////////////
    for(int iW = 0; iW < nSamples; iW++){
        const STATE currentW = (iW+1) * wSample;
        TPZAutoPointer<TPZStructMatrix> structMatrixPtr = an.StructMatrix();
        structMatrixPtr->SetMaterialIds(matsNonPml);
        //assemble load vector
        an.AssembleResidual();
        //assemble pml
        structMatrixPtr->SetMaterialIds(matsPml);
        for(auto itMap : matsPml){
            TPZMatAcousticsPml *mat = dynamic_cast<TPZMatAcousticsPml *>(cmesh->FindMaterial(itMap));
            if(mat!=nullptr){
                mat->SetW(currentW);
            }else{
                DebugStop();
            }
            an.Assemble();
            //sum all matrices
            //solve system
            //get solution
            //inverse transform
        }
    }
    //loads real solution
    an.LoadSolution();

//    boost::posix_time::ptime t1_a =
//            boost::posix_time::microsec_clock::local_time();
//    an.Assemble();
//    boost::posix_time::ptime t2_a =
//            boost::posix_time::microsec_clock::local_time();
//    boost::posix_time::ptime t1_s =
//            boost::posix_time::microsec_clock::local_time();
//    an.Solve();
//    boost::posix_time::ptime t2_s =
//            boost::posix_time::microsec_clock::local_time();
//
//
//    std::cout<<"Assembly duration: "<<t2_a-t1_a
//             <<"  Solver duration: "<<t2_s-t1_s<<std::endl;
//    an.LoadSolution();
//    if (l2error) {
//        std::cout<<"Calculating errors..."<<std::endl;
//        TPZVec<REAL> errorVec(1,0);
//        an.SetExact(&exactSol);
//        errorVec.Resize(2, 0.);
//        an.SetThreadsForError(nThreads);
//        an.PostProcessError(errorVec);
//        std::string fileName = prefix + "error_nel";
//        fileName.append(std::to_string(nDiv * nDiv * 2));
//        fileName.append("_p");
//        fileName.append(std::to_string(pOrder));
//        fileName.append(".csv");
//        std::ofstream errorFile(fileName.c_str());
//        errorFile << neq << "," << errorVec[0] << "," << errorVec[1] << "," << errorVec[2]
//                  << std::endl;
//        errorFile.close();
//        std::cout<<" Done!"<<std::endl;
//    }

    if (genVTK) {
        std::cout<<"Post processing... ";
        TPZStack<std::string> scalnames, vecnames;
        vecnames.Push("E");
        std::string plotfile = prefix+"sol";
        plotfile.append(std::to_string(cmesh->NElements()));
        plotfile.append(".vtk");

        an.DefineGraphMesh(2, scalnames, vecnames,
                           plotfile);  // define malha grafica
        int postProcessResolution = postprocessRes; // define resolucao do pos processamento
        an.PostProcess(postProcessResolution);
        std::cout<<" Done!"<<std::endl;
    }
    gmesh->SetReference(nullptr);
    cmesh->SetReference(nullptr);
    delete cmesh;
    delete gmesh;
}

void
CreateGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const bool &print,
            const REAL &elSize, const REAL &length, const REAL &height, const REAL &pmlLength, const std::string &prefix) {

    std::ostringstream str_elSize;
    str_elSize << std::setprecision(20) << elSize;
    std::ostringstream str_length;
    str_length<< std::setprecision(20) << length;
    std::ostringstream str_height;
    str_height<< std::setprecision(20) << height;
    std::ostringstream str_pmlLength;
    str_pmlLength<< std::setprecision(20) << pmlLength;

    std::string command = "gmsh " + mshFileName + " -2 -match ";
    command += " -v 3 ";
    command += " -setnumber length "+str_length.str();
    command += " -setnumber height "+str_height.str();
    command += " -setnumber el_size "+str_elSize.str();
    command += " -setnumber pml_length "+str_pmlLength.str();
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
    auto matIds = meshReader.fMaterialDataVec.fMatID;
    matIdVec.Resize(matIds.size());
    int i = 0;
    for(auto id = matIds.begin(); id != matIds.end(); id++,i++ )    {
        matIdVec[i] = *id;
    }
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
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal) {
    TPZManVector<long, 1000> allConnects;
    std::set<long> boundConnects;

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
            std::set<long> boundConnectsEl;
            cel->BuildConnectList(boundConnectsEl);

            for (std::set<long>::iterator iT = boundConnectsEl.begin();
                 iT != boundConnectsEl.end(); iT++) {
                const long val = *iT;
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
CreateCMesh(TPZCompMesh *&cmesh, TPZGeoMesh *gmesh, int pOrder, void (&loadVec)(const TPZVec<REAL> &, TPZVec<STATE> &),
            const REAL &alphaPml, const std::string &prefix, const bool &print,
            const TPZVec<int> &matIdVec, const REAL &rho, const REAL &velocity) {
    const int dim = 2;   // dimensao do problema
    const int matId = 1; // define id para um material(formulacao fraca)
    const int bc0 = -1;  // define id para um material(cond contorno dirichlet)
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
    TPZMatAcousticsH1 *matAcoustics = NULL;
    matAcoustics = new TPZMatAcousticsH1(matIdVec[0], rho, velocity);
    matAcoustics->SetForcingFunction(loadVec, 4);
    cmesh->InsertMaterialObject(matAcoustics);

    TPZMatAcousticsPml *matAcousticsPML = NULL;
    REAL pmlBegin, pmlLength;
    {
        REAL xMax =-1e20,xMin = 1e20;
        TPZGeoMesh * gmesh = cmesh->Reference();
        for (int iel = 0; iel < gmesh->NElements(); ++iel) {
            TPZGeoEl *geo = gmesh->Element(iel);
            if (geo->MaterialId() == matIdVec[1]) {
                for (int iNode = 0; iNode < geo->NCornerNodes(); ++iNode) {
                    TPZManVector<REAL, 3> co(3);
                    geo->Node(iNode).GetCoordinates(co);
                    const REAL &xP = co[0];
                    const REAL &yP = co[1];
                    if (xP > xMax) {
                        xMax = xP;
                    }
                    if (xP < xMin) {
                        xMin = xP;
                    }
                }
            }
        }
        pmlBegin = xMin;
        pmlLength = xMax - xMin;
    }
    //matAcousticsPML = new TPZMatAcousticsPml(matIdVec[1], rho, velocity, pmlBegin, pmlLength, alphaPml);
    matAcousticsPML =
            new TPZMatAcousticsPml(matIdVec[1],*matAcoustics,true, pmlBegin, false, -1, alphaPml,pmlLength);
    cmesh->InsertMaterialObject(matAcousticsPML);
    {
        REAL xMax =-1e20,xMin = 1e20;
        TPZGeoMesh * gmesh = cmesh->Reference();
        for (int iel = 0; iel < gmesh->NElements(); ++iel) {
            TPZGeoEl *geo = gmesh->Element(iel);
            if (geo->MaterialId() == matIdVec[2]) {
                for (int iNode = 0; iNode < geo->NCornerNodes(); ++iNode) {
                    TPZManVector<REAL, 3> co(3);
                    geo->Node(iNode).GetCoordinates(co);
                    const REAL &xP = co[0];
                    const REAL &yP = co[1];
                    if (xP > xMax) {
                        xMax = xP;
                    }
                    if (xP < xMin) {
                        xMin = xP;
                    }
                }
            }
        }
        pmlBegin = xMax;
        pmlLength = xMax - xMin;
    }
    matAcousticsPML =
            new TPZMatAcousticsPml(matIdVec[1],*matAcoustics,true, pmlBegin, false, -1, alphaPml,pmlLength);
    cmesh->InsertMaterialObject(matAcousticsPML);

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
        std::string fileName(prefix+ "cmesh");
        fileName.append(".txt");
        std::ofstream fileHCurl(fileName.c_str());
        cmesh->Print(fileHCurl);
    }
}
