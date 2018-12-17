#include "TPZAcousticGeoMesher.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeopoint.h"
#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "pzcheckgeom.h"
#include "pzgeoelrefless.h"

#include <array>


void TPZAcousticGeoMesher::GetMaterialProperties(std::map<std::string, int> &materialNames, std::map<int, REAL> &rhoMap,
                                                 std::map<int, REAL> &velocityMap) {
    std::map<std::string,std::pair<REAL,REAL>> dataBase;//rho,velocity
    dataBase.insert(std::make_pair<>("water",std::make_pair<REAL,REAL>(1000,1500)));
    dataBase.insert(std::make_pair<>("steel",std::make_pair<REAL,REAL>(7850,5960)));
    dataBase.insert(std::make_pair<>("cement",std::make_pair<REAL,REAL>(3150,3500)));
    dataBase.insert(std::make_pair<>("rock",std::make_pair<REAL,REAL>(2750,5950)));

    for(auto &itMap : materialNames){
        auto possibleMat = dataBase.find(itMap.first);
        if(possibleMat == dataBase.end()){
            PZError << "You are trying to create a nonexistent material called "<< itMap.first <<std::endl;
            PZError << "Execution will abort now."<<std::endl;
            DebugStop();
        }
        rhoMap[itMap.second] = possibleMat->second.first;
        velocityMap[itMap.second] = possibleMat->second.second;
    }
}

void TPZAcousticGeoMesher::ReadMeshMaterials(const std::string &fileName, const std::string &prefix,
                                             std::map<int, REAL> &velocityMap, std::map<int, REAL> &rhoMap) {
    std::string command = "gmsh " + fileName + " -0 -v 0 -format msh2";
    command += " -o " + prefix + "wellMesh.msh";
    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");

    TPZGmshReader meshReader;
    TPZGeoMesh *gmesh = meshReader.GeometricGmshMesh(prefix+"wellMesh.msh", nullptr, false);
    std::map<std::string,int> materialNames = meshReader.fPZMaterialId[2];
    GetMaterialProperties(materialNames, rhoMap, velocityMap);
    delete gmesh;
}
void TPZAcousticGeoMesher::ReadGmshMesh(TPZGeoMesh * &gmesh, const std::string &meshFileName,
        const std::string &prefix, const TPZVec<REAL> &elSizes, TPZVec<std::map<int,int>> &translatedMatIds){
    std::string command = "gmsh " + meshFileName + " -2 -match -format msh2";
    command += " -v 3 ";

    for(int i = 0; i < elSizes.size(); i++){
        std::ostringstream str_elSize;
        str_elSize << std::setprecision(16) << elSizes[i];
        command += " -setnumber el_size_"+std::to_string(i+1)+" "+str_elSize.str();
    }
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
    translatedMatIds = meshReader.fMatIdTranslate;
}
void TPZAcousticGeoMesher::CreateGMesh(TPZGeoMesh *&gmesh, const std::string meshFileName, const std::string &prefix,
                                       REAL nElemPerLambdaTimesOmega, TPZVec<int> &matIdVec, TPZVec <REAL> &elSizes,
                                       std::map<int, REAL> &velocityMap, std::map<int, REAL> &rhoMap) {
    ReadMeshMaterials(meshFileName, prefix, velocityMap, rhoMap);
    int i =0;
    for(auto iVelocity : velocityMap){
        elSizes[i] = 2 *M_PI*iVelocity.second / nElemPerLambdaTimesOmega;
        i++;
    }
    TPZManVector<std::map<int,int>,5> translatedMatIds;
    ReadGmshMesh(gmesh, meshFileName, prefix, elSizes, translatedMatIds);
#ifdef PZDEBUG
    TPZCheckGeom * Geometrytest = new TPZCheckGeom(gmesh);
    int isBadMeshQ = Geometrytest->PerformCheck();

    if (isBadMeshQ) {
        DebugStop();
    }
#endif
    const int n1dMat =translatedMatIds[1].size();
    const int n2dMat =translatedMatIds[2].size();
    matIdVec.Resize(n2dMat+n1dMat);
    int imat = 0;
    ///at this point, the materials in the .geo must be declared in a crescent order
    for (auto& kv :translatedMatIds[2]) {
        matIdVec[imat] = kv.first;
        imat++;
    }
    for (auto& kv :translatedMatIds[1]) {
        matIdVec[imat] = kv.first;
        imat++;
    }
    return;
}

void TPZAcousticGeoMesher::CreateSourceNode(TPZGeoMesh *&gmesh, const int &matIdSource, const REAL &sourcePosX,
                                            const REAL &sourcePosY) {
    int64_t sourceNodeIndex = -1;
    REAL minDist = 1e12;
    //////
    for(int iEl = 0; iEl < gmesh->NElements(); iEl++){
        const TPZGeoEl & currentEl = *(gmesh->ElementVec()[iEl]);
        const int nNodes = currentEl.NCornerNodes();
        for(int iNode = 0; iNode < nNodes; iNode++){
            const TPZGeoNode & currentNode = currentEl.Node(iNode);
            const REAL xNode = currentNode.Coord(0);
            const REAL yNode = currentNode.Coord(1);
            const REAL currentDist =
                    (sourcePosX - xNode)*(sourcePosX - xNode) + (sourcePosY - yNode)*(sourcePosY - yNode);
            if(currentDist < minDist){
                minDist = currentDist;
                sourceNodeIndex = currentNode.Id();
            }
        }
    }
    TPZVec<int64_t> nodeIdVec(1,sourceNodeIndex);
    TPZGeoElRefLess<pzgeom::TPZGeoPoint > *zeroDEl =
            new TPZGeoElRefLess<pzgeom::TPZGeoPoint >(nodeIdVec, matIdSource,
                                                      *gmesh);
    gmesh->BuildConnectivity();
}

void TPZAcousticGeoMesher::PrintMesh(TPZGeoMesh *&gmesh, const std::string &fileName, const std::string &prefix,
                                     bool printVTK, bool printTxt) {
    const std::string meshFileName = prefix + fileName;
    if(printVTK){
        const std::string vtkFile = meshFileName + ".vtk";
        std::ofstream outVTK(vtkFile.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
    if(printTxt){
        const std::string txtFile = meshFileName + ".vtk";
        std::ofstream outTXT(txtFile.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
}
