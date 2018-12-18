#ifndef TPZACOUSTICGEOMESHER
#define TPZACOUSTICGEOMESHER


#include <map>
#include "pzvec.h"
class TPZGeoMesh;

/**
 * A class that interfaces with gmsh executable  in the context of acoustic simulations. It also handles the
 * generation of a NeoPZ geometric mesh.
 *
 */
class TPZAcousticGeoMesher{
public:

    /**
     * This constructor will not be generated
     */
    TPZAcousticGeoMesher() = delete;

    /**
     * This initializes the TPZAcousticGeoMesher structure
     * @param meshFileName name of the .geo file (in)
     * @param prefix relative path to the .geo file (in)
     */
    TPZAcousticGeoMesher(const std::string meshFileName, const std::string &prefix);
    /**
     * This method generates a TPZGeoMesh from a .geo file. It returns all the relevant information regarding the
     * created mesh.
     * @param nElemPerLambdaTimesOmega how many elements are created for each wavelength (times ), i.e.,
     * elSize = 2 pi v /(nElemPerLambdaTimesOmega) = 2 pi v/(nelem 2pi f) = lambda / nelem (in)
     * @param matIdVec vector with the material ids in the order they were declared (as physical surfaces)
     * in the .geo file (out)
     */
    void CreateGMesh(REAL nElemPerLambdaTimesOmega, TPZVec<int> &matIdVec);
    /**
     * Creates a 0D element (node) representing the source. The location is approximated, since the source will be
     * placed at the nearest node in the mesh, so there is a maximum error of elSize/2
     * @param gmesh the mesh in which the source will be inserted (in/out)
     * @param matIdSource material id associated with the source (in)
     * @param sourcePosX x-coordinate for the source (in)
     * @param sourcePosY y-coordinate for the source (in)
     */
    void CreateSourceNode(const int &matIdSource, const REAL &sourcePosX, const REAL &sourcePosY);

    /**
     * Prints the geometric mesh in both .vtk and .txt formats.
     * @param gmesh mesh to be printed (in)
     * @param fileName where to print the mesh (extensions will be added to fileName)
     * @param prefix directory to output the printed mesh (in)
     * @param printVTK whether to print the mesh in .vtk format (in)
     * @param printTxt whether to print the mesh in .txt format (in)
     */
    void PrintMesh(const std::string &fileName, const std::string &prefix, bool printVTK, bool printTxt);
    /**
     *
     * @return Materials' velocities indexed by material id
     */
    std::map<int,REAL> GetVelocityMap() const { return fVelocityMap;}
    /**
     *
     * @return Materials' densities indexed by material id
     */
     std::map<int,REAL> GetDensityMap() const { return fRhoMap;}

     /**
      *
      * @return Materials' element sizes indexed by material id
      */
     std::map<int,REAL> GetElSizes() const;

     TPZGeoMesh *GetMesh() const{return fGmesh; }
protected:
    void ReadMeshMaterials();
    /**
     * Generates a .msh (msh version 2) from a .geo file (calling the gmsh executable) and, finally, creates
     * a TPZGeoMesh from the .msh file.
     * @param gmesh mesh to be created (out)
     * @param meshFileName name of the .geo file (in)
     * @param prefix directory containing the .geo file (in)
     * @param elSizes vector with element sizes for each material. the .geo file must have variables declared as
     * DefineConstant[el_size_i = 0.0125]; where i = 1,2,3,... (in the same order that the physical materials
     * are declared)
     */
    void ReadGmshMesh(TPZVec<std::map<int,int>> &translatedMatIds);

    /**
     * This method will fill two maps with the properties (sound speed and density) of the materials
     * present in the mesh.
     * @param materialNames material names indexed by their material ids (in)
     * @param rhoMap material densities indexed by their material ids (out)
     * @param velocityMap material velocities indexed by their material ids (out)
     */
    void GetMaterialProperties(std::map<std::string,int> &materialNames,std::map<int,REAL> &rhoMap,
                                      std::map<int,REAL> &velocityMap);

    /**
     * This points to the actual mesh
     */
    TPZGeoMesh *fGmesh;
    /**
     * Name of the .geo file on which the mesh will be based
     */
    const std::string fMeshFileName;
    /**
     * Relative path to the .geo file indicated by fMeshFileName
     */
    const std::string fPrefix;

    /**
     * Material velocities indexed by their material ids
     */
    std::map<int, REAL> fVelocityMap;
    /**
     * Material densities indexed by their material ids
     */
     std::map<int, REAL> fRhoMap;
     /**
      * Vector containing the elSize for each material
      */
     TPZVec<REAL> fElSizes;
};

#endif //TPZACOUSTICGEOMESHER