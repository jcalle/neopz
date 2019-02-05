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

void TPZAcousticGeoMesher::ReadMeshMaterials() {
    std::string command = "gmsh " + fMeshFileName+ " -0 -v 0 -format msh2";
    command += " -o " + fPrefix + "wellMesh.msh";
    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");

    TPZGmshReader meshReader;
    TPZGeoMesh *gmesh = meshReader.GeometricGmshMesh(fPrefix+"wellMesh.msh", nullptr, false);
    std::map<std::string,int> materialNames = meshReader.fPZMaterialId[2];
    GetMaterialProperties(materialNames, fRhoMap, fVelocityMap);
    delete gmesh;
}

void TPZAcousticGeoMesher::ReadGmshMesh(TPZVec<std::map<int,int>> &translatedMatIds){
    std::string command = "gmsh " + fMeshFileName + " -2 -match -format msh2";
    command += " -v 3 ";

    for(int i = 0; i < fElSizes.size(); i++){
        std::ostringstream str_elSize;
        str_elSize << std::setprecision(16) << fElSizes[i];
        command += " -setnumber el_size_"+std::to_string(i+1)+" "+str_elSize.str();
    }
    command += " -o " + fPrefix + "wellMesh.msh";
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
    fGmesh = meshReader.GeometricGmshMesh(fPrefix+"wellMesh.msh");
    translatedMatIds = meshReader.fMatIdTranslate;
}
void TPZAcousticGeoMesher::CreateGMesh(REAL nElemPerLambdaTimesOmega) {
    if(fGmesh!=nullptr){
        std::cout<<"You have requested to overwrite the mesh, so the existing mesh will be deleted first"<<std::endl;
        delete fGmesh;
    }
    int i =0;
    fElSizes.Resize( fVelocityMap.size() );
    for(auto iVelocity : fVelocityMap){
        fElSizes[i] = 2 *M_PI*iVelocity.second / nElemPerLambdaTimesOmega;
        i++;
    }
    TPZManVector<std::map<int,int>,5> translatedMatIds;
    ReadGmshMesh(translatedMatIds);
#ifdef PZDEBUG
    TPZCheckGeom * Geometrytest = new TPZCheckGeom(fGmesh);
    int isBadMeshQ = Geometrytest->PerformCheck();

    if (isBadMeshQ) {
        DebugStop();
    }
#endif
    const int n1dMat =translatedMatIds[1].size();
    const int n2dMat =translatedMatIds[2].size();
    fMatIdVec.Resize(n2dMat+n1dMat);
    int imat = 0;
    ///at this point, the materials in the .geo must be declared in a crescent order
    for (auto& kv :translatedMatIds[2]) {
        fMatIdVec[imat] = kv.first;
        imat++;
    }
    for (auto& kv :translatedMatIds[1]) {
        fMatIdVec[imat] = kv.first;
        imat++;
    }
    return;
}

void TPZAcousticGeoMesher::CreateSourceNode(const REAL &sourcePosX,
                                            const REAL &sourcePosY) {
#ifdef PZDEBUG
    if(fGmesh == nullptr){
        PZError<<"You are trying to create a source node in an empty mesh. Aborting..."<<std::endl;
        DebugStop();
    }
#endif

    for (int iMat = 0; iMat< fMatIdVec.size(); iMat++){
        fMatIdSource += fMatIdVec[iMat];
    }
    TPZVec<int> matIdVecCp(fMatIdVec);
    fMatIdVec.resize(fMatIdVec.size()+1);
    fMatIdVec[0] = fMatIdSource;
    for (int iMat = 0; iMat< matIdVecCp.size(); iMat++){
        fMatIdVec[iMat+1] = matIdVecCp[iMat];
    }


    int64_t sourceNodeIndex = -1;
    REAL minDist = 1e12;
    //////
    for(int iEl = 0; iEl < fGmesh->NElements(); iEl++){
        const TPZGeoEl & currentEl = *(fGmesh->ElementVec()[iEl]);
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
            new TPZGeoElRefLess<pzgeom::TPZGeoPoint >(nodeIdVec, fMatIdSource,
                                                      *fGmesh);
    fGmesh->BuildConnectivity();
}

void TPZAcousticGeoMesher::PrintMesh(const std::string &fileName, const std::string &prefix,
                                     bool printVTK, bool printTxt) {
    const std::string meshFileName = prefix + fileName;
    if(printVTK){
        const std::string vtkFile = meshFileName + ".vtk";
        std::ofstream outVTK(vtkFile.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(fGmesh, outVTK, true);
        outVTK.close();
    }
    if(printTxt){
        const std::string txtFile = meshFileName + ".vtk";
        std::ofstream outTXT(txtFile.c_str());
        fGmesh->Print(outTXT);
        outTXT.close();
    }
}

TPZAcousticGeoMesher::TPZAcousticGeoMesher(const std::string meshFileName, const std::string &prefix) :
fMeshFileName(meshFileName) ,fPrefix(prefix), fMatIdSource(-1) {
    fGmesh = nullptr;
    fMatIdVec.Resize(0);
    ReadMeshMaterials();
}

std::map<int, REAL> TPZAcousticGeoMesher::GetElSizes() const {
    std::map<int,REAL> elSizeMap;
    int i = 0;
    for(auto itMap : fVelocityMap){
        elSizeMap[itMap.first] = fElSizes[i];
        i++;
    }
    return elSizeMap;
}
