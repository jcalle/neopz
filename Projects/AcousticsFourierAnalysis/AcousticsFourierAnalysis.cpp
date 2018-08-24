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
#include <pzl2projection.h>
#include <TPZSkylineNSymStructMatrix.h>
#include <pzskylnsymmat.h>
#include <pzgeopoint.h>
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

void FilterBoundaryEquations(TPZCompMesh *cmeshHCurl,
                             TPZVec<int64_t> &activeEquations, int &neq,
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
    int pOrder = 2; //PARAMS
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

void RunSimulation(const int &nDiv, const int &pOrder, const std::string &prefix, const REAL &wZero) {
    // PARAMETROS FISICOS DO PROBLEMA
    const int nThreads = 8; //PARAMS
    const bool l2error = true; //PARAMS
    const bool genVTK = true; //PARAMS
    const bool printG = true;//PARAMS
    const bool printC = true;//PARAMS
    const int postprocessRes = 0;//PARAMS
    REAL alphaPML = 10;
    REAL rho = 1.3;
    REAL velocity = 340;
    REAL peakTime = 1./100;
    REAL amplitude = 10;
    REAL totalTime = 5 * peakTime;
    const int64_t nTimeSteps = 100;
    REAL elSize = 2 *M_PI*velocity / (12 *wZero),length = 60,height = 5,pmlLength = 20;



    std::cout<<"Creating gmesh... ";
    boost::posix_time::ptime t1_g =
            boost::posix_time::microsec_clock::local_time();
    TPZGeoMesh *gmesh = NULL;
    std::string meshFileName("wellMesh.geo");
    TPZVec<int> matIdVec;
    CreateGMesh(gmesh, meshFileName, matIdVec, printG, elSize, length,
                height, pmlLength, prefix);

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

    {
        const REAL sourcePosX = length/2.;
        const REAL sourcePosY = height/2.;
        REAL minDist = 1e6;
        int64_t nodeIndex = -1;
        for(int iNode = 0; iNode < gmesh->NNodes(); iNode++){
            const TPZGeoNode & currentNode = gmesh->NodeVec()[iNode];
            const REAL xNode = currentNode.Coord(0);
            const REAL yNode = currentNode.Coord(1);
            const REAL currentDist =
                    (sourcePosX - xNode)*(sourcePosX - xNode) + (sourcePosY - yNode)*(sourcePosY - yNode);
            if(currentDist < minDist){
                minDist = currentDist;
                nodeIndex = currentNode.Id();
            }
        }
        TPZVec<int64_t> nodeIdVec(1,nodeIndex);
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
    CreateCMesh(cmesh, gmesh, pOrder, loadVec, alphaPML, prefix, printC, matIdVec, rho, velocity);
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
    TPZSkylineNSymStructMatrix structMatrix(cmesh);
#else
    TPZSpStructMatrix structMatrix(cmesh);
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
    ////////////////////////////////////////////////////////////////////////
    const REAL wMax = 3*wZero;
    const int nSamples = 350;
    const REAL wSample = wMax/nSamples;

    ////////////////////////////////////////////////////////////////////////
    std::set<int> matsPml;
    std::set<int> matsNonPml;
    std::set<int> matSource;
    for(auto itMap : cmesh->MaterialVec()) {
        TPZMatAcousticsPml *matPml = dynamic_cast<TPZMatAcousticsPml *>(itMap.second);
        TPZMatAcousticsH1 *matH1 = dynamic_cast<TPZMatAcousticsH1 *>(itMap.second);
        if (matH1 != nullptr && matPml == nullptr) {
            if(matH1->Id() == matIdVec[0]){//TODO: fazer isso de maneira menos tosca. queremos pular o elemento da fonte
                matSource.insert(matH1->Id());
                continue;
            }
            matsNonPml.insert(matH1->Id());
        } else if (matPml != nullptr) {
            matsPml.insert(matPml->Id());
        }
    }


#ifdef USING_SKYLINE
    TPZSkylNSymMatrix<STATE> matM;
    TPZSkylNSymMatrix<STATE> matK;
#else
    TPZFYsmpMatrix<STATE> matM;
    TPZFYsmpMatrix<STATE> matK;
#endif


    {
        TPZAutoPointer<TPZStructMatrix> structMatrixPtr = an.StructMatrix();
        structMatrixPtr->SetMaterialIds(matsNonPml);
        //set structMatrix M
        for (std::set<int>::iterator it=matsNonPml.begin(); it!=matsNonPml.end(); ++it){
            TPZMatAcousticsH1 *matH1 = dynamic_cast<TPZMatAcousticsH1 *>(cmesh->MaterialVec()[*it]);
            matH1->SetAssemblingMatrix(TPZMatAcousticsH1::M);
        }
        an.Assemble();
        //get structMatrix M

#ifdef USING_SKYLINE
        auto mat = dynamic_cast<TPZSkylNSymMatrix<STATE> *>( an.Solver().Matrix().operator->() );
#else
        auto mat = dynamic_cast<TPZFYsmpMatrix<STATE> *>( an.Solver().Matrix().operator->() );
#endif
        if(!mat){
            DebugStop();
        }
        matM = *mat;
        //set structMatrix K
        for (std::set<int>::iterator it=matsNonPml.begin(); it!=matsNonPml.end(); ++it){
            TPZMatAcousticsH1 *matH1 = dynamic_cast<TPZMatAcousticsH1 *>(cmesh->MaterialVec()[*it]);
            matH1->SetAssemblingMatrix(TPZMatAcousticsH1::K);
        }
        an.Assemble();
        //get structMatrix K
#ifdef USING_SKYLINE
        mat = dynamic_cast<TPZSkylNSymMatrix<STATE> *>( an.Solver().Matrix().operator->() );
#else
        mat = dynamic_cast<TPZFYsmpMatrix<STATE> *>( an.Solver().Matrix().operator->() );
#endif
        if(!mat){
            DebugStop();
        }
        matK = *mat;
    }


    ////////////////////////////////////////////////////////////////////////
    TPZFMatrix<STATE> timeDomainSolution(neq,nTimeSteps,0.);
    TPZFMatrix<STATE> rhsFake(cmesh->NEquations(),1);
#ifdef USING_SKYLINE
    TPZAutoPointer<TPZMatrix<STATE>> matFinal(new TPZSkylNSymMatrix<STATE>(matK));
#else
    TPZAutoPointer<TPZMatrix<STATE>> matFinal(new TPZFYsmpMatrix<STATE>(matK));
#endif
    for(int iW = 0; iW < nSamples - 1; iW++){
        const STATE currentW = (iW+1) * wSample;
        TPZAutoPointer<TPZStructMatrix> structMatrixPtr = an.StructMatrix();
        //assemble pml
        structMatrixPtr->SetMaterialIds(matsPml);
        for(auto itMap : matsPml) {
            TPZMatAcousticsPml *mat = dynamic_cast<TPZMatAcousticsPml *>(cmesh->FindMaterial(itMap));
            if (mat != nullptr) {
                mat->SetW(currentW);
            } else {
                DebugStop();
            }
        }

        //-w^2 M + K
        boost::posix_time::ptime t1_sum =
          boost::posix_time::microsec_clock::local_time();



#ifdef USING_SKYLINE
        for(int iCol = 0; iCol < neq; iCol++){
            TPZSkylNSymMatrix<STATE> * matFinalDummy = dynamic_cast<TPZSkylNSymMatrix<STATE> *>(matFinal.operator->());
            const int nRows = matFinalDummy->SkyHeight(iCol);
            for(int iRow = iCol; iRow >= iCol - nRows; iRow --){
                matFinal->PutVal(iRow,iCol, -1.*currentW*currentW*matM.GetVal(iRow,iCol) + matK.GetVal(iRow,iCol));
                matFinal->PutVal(iCol,iRow, -1.*currentW*currentW*matM.GetVal(iCol,iRow) + matK.GetVal(iCol,iRow));
            }
        }
#else
        TPZVec<int64_t> ia;
        TPZVec<int64_t> ja;
        TPZVec<STATE> aM, aK, aFinal;
        matM.GetData(ia,ja,aM);
        matK.GetData(ia,ja,aK);
        aFinal.Resize(aK.size());

        TPZFYsmpMatrix<STATE> * matFinalDummy = dynamic_cast<TPZFYsmpMatrix<STATE> *>(matFinal.operator->());
        for(int i = 0; i < aFinal.size(); i++){
            aFinal[i] = -1.*currentW*currentW*aM[i] + aK[i];
        }
        matFinalDummy->SetData(ia,ja,aFinal);
#endif


        boost::posix_time::ptime t2_sum =
          boost::posix_time::microsec_clock::local_time();
        std::cout<<"Summing matrices took "<<t2_sum-t1_sum<<std::endl;
        //pml elements
        structMatrixPtr->Assemble(matFinal,rhsFake,nullptr);


        structMatrixPtr->SetMaterialIds(matSource);

        TPZMatAcousticsH1 *matSourcePtr = dynamic_cast<TPZMatAcousticsH1 *>(cmesh->MaterialVec()[matIdSource]);
        if(matSourcePtr == nullptr){
            DebugStop();
        }
//        auto currentSource = [currentW, amplitude, peakTime, wZero](const TPZVec<REAL> &coord, TPZVec<STATE> &result) {
//            result.Resize(1, 0.);
//            result[0] += -2. * sqrt(-1) * amplitude * currentW;
//            result[0] *= exp(0.5 + sqrt(-1) * currentW * peakTime - (currentW/wZero)*(currentW/wZero));
//            result[0] *= sqrt(2 * M_PI)/( wZero * wZero);
//        };

        STATE currentSource = 0;
        currentSource += -2. * sqrt((STATE)-1) * amplitude * currentW;
        currentSource *= exp(0.5 + sqrt((STATE)-1) * currentW * peakTime - (currentW/wZero)*(currentW/wZero));
        currentSource *= sqrt(2 * M_PI)/( wZero * wZero);
        matSourcePtr->SetSourceFunc(currentSource);
        //assemble load vector
        an.AssembleResidual();
        //solve system
        an.Solver().SetMatrix(matFinal);

        std::cout<<"solver step "<<iW+1<<" out of "<<nSamples - 1;
        boost::posix_time::ptime t1_solve =
          boost::posix_time::microsec_clock::local_time();
        std::cout<<std::flush;
        an.Solve();
        boost::posix_time::ptime t2_solve =
          boost::posix_time::microsec_clock::local_time();
        std::cout<<"solver step "<<iW+1<<" out of "<<nSamples - 1<<" took "<<t2_solve-t1_solve<<std::endl;
        //get solution
        TPZFMatrix<STATE> &currentSol = an.Solution();//p omega
        //inverse transform

        REAL timeStepSize = totalTime/nTimeSteps;
        for(int iPt = 0; iPt < neq; iPt++){
            for(int iTime = 0; iTime < nTimeSteps; iTime++){
                timeDomainSolution(iPt,iTime) +=
                        std::real(currentSol(iPt,0)) * wSample * std::cos((REAL)iTime * timeStepSize * currentW)/M_PI;
            }
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
            std::cout<<"\rtime: "<<iTime+1<<" out of "<<nTimeSteps<<std::endl;
            std::cout<<std::flush;
            for(int iPt = 0; iPt < neq; iPt++){
                currentSol(iPt,0) = timeDomainSolution(iPt,iTime);
            }
            if(filter){
                structMatrix.EquationFilter().Scatter(currentSol, scatteredSol);
                an.LoadSolution(scatteredSol);
            }else{
                an.LoadSolution(currentSol);
            }
            an.PostProcess(postProcessResolution);
        }
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
    auto matIds = meshReader.fMaterialDataVec;
    matIdVec.Resize(4);
    //fPZMaterialId[dimension][name]
    matIdVec[0] = meshReader.fPZMaterialId[2]["water"];
    matIdVec[1] = meshReader.fPZMaterialId[2]["pmlLeft"];
    matIdVec[2] = meshReader.fPZMaterialId[2]["pmlRight"];
    matIdVec[3] = meshReader.fPZMaterialId[1]["bound"];

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
CreateCMesh(TPZCompMesh *&cmesh, TPZGeoMesh *gmesh, int pOrder, void (&loadVec)(const TPZVec<REAL> &, TPZVec<STATE> &),
            const REAL &alphaPml, const std::string &prefix, const bool &print,
            const TPZVec<int> &matIdVec, const REAL &rho, const REAL &velocity) {
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
    TPZMatAcousticsH1 *matAcoustics = NULL;
    matAcoustics = new TPZMatAcousticsH1(matIdVec[0], rho, velocity);
    matAcoustics->SetForcingFunction(loadVec, 4);
    cmesh->InsertMaterialObject(matAcoustics);
    matAcoustics = new TPZMatAcousticsH1(matIdVec[1], rho, velocity);
//    matAcoustics->SetForcingFunction(loadVec, 4);
    cmesh->InsertMaterialObject(matAcoustics);
//    auto exactSol = [](const TPZVec<REAL> &coord, TPZVec<STATE> &result, TPZFMatrix<STATE> &grad) {
//        result.Resize(1, 0.);
//        result[0] = M_PI * cos((1/10.) * 2 * M_PI * coord[0]) * sin((1/10.) * 2 * M_PI  * coord[1]);
//
//        grad.Resize(2, 1);
//        grad(0, 0) = -1 * M_PI * M_PI * sin(M_PI * coord[0]) * sin(M_PI * coord[1]);
//        grad(1, 0) = M_PI * M_PI * cos(M_PI * coord[0]) * cos(M_PI * coord[1]);
//    };
//    matAcoustics->SetExactSol(exactSol);

    TPZMatAcousticsPml *matAcousticsPML = NULL;
    REAL pmlBegin, pmlLength;
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
        pmlBegin = xMin;
        pmlLength = xMax - xMin;
    }
    //matAcousticsPML = new TPZMatAcousticsPml(matIdVec[1], rho, velocity, pmlBegin, pmlLength, alphaPml);
    matAcousticsPML =
            new TPZMatAcousticsPml(matIdVec[2],*matAcoustics,true, pmlBegin, false, -1, alphaPml,pmlLength);
//    matAcousticsPML->SetExactSol(exactSol);
    cmesh->InsertMaterialObject(matAcousticsPML);
    {
        REAL xMax =-1e20,xMin = 1e20;
        TPZGeoMesh * gmesh = cmesh->Reference();
        for (int iel = 0; iel < gmesh->NElements(); ++iel) {
            TPZGeoEl *geo = gmesh->Element(iel);
            if (geo->MaterialId() == matIdVec[3]) {
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
            new TPZMatAcousticsPml(matIdVec[3],*matAcoustics,true, pmlBegin, false, -1, alphaPml,pmlLength);
//    matAcousticsPML->SetExactSol(exactSol);
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
