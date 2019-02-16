//
//  TPZGmshReader.h
//  PZ
//
//  Created by Omar on 2/7/16.
//
//

#ifndef TPZGmshReader_h
#define TPZGmshReader_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "pzgmesh.h"


class TPZGeoMesh;


struct MaterialDataS {
    
    TPZStack<int> fMatID;
    TPZStack<std::pair<int ,std::string> >  fMaterial;
    
    MaterialDataS() : fMatID(), fMaterial(){
        
    }
    
    MaterialDataS(int num) : fMatID(), fMaterial(){
        
    }
    
    MaterialDataS(const MaterialDataS &copy) : fMatID(copy.fMatID),
    fMaterial(copy.fMaterial) {
    }
    
    MaterialDataS &operator=(const MaterialDataS &copy){
        fMatID = copy.fMatID;
        fMaterial = copy.fMaterial;
        return *this;
    }
    
};


/**
 * @brief Implement the interface between TPZGeoMesh and the files produced by Gmsh (version 3.0 or 4.0 ) in msh format.
 * @since January 16, 2017
 */

/** What is Gmsh ? Take a look on http://gmsh.info/
 * Gmsh is a free 3D finite element grid generator with a build-in CAD engine and post-processor. Its design goal is to provide a fast, light and user-friendly
 * meshing tool with parametric input and advanced visualization capabilities. Gmsh is built around four modules: geometry, mesh, solver and post-processing.
 * The specification of any input to these modules is done either interactively using the graphical user interface or in ASCII text files using Gmsh's own
 * scripting language.
 */

/** Note about the implementation for file format 4
 * The mandatory sections are considered MeshFormat, Entities, Nodes and Elements.
 * The optional section PhysicalName is considered the others (PartitionedEntities,Periodic,GhostElements,NodeData,ElementData,ElementNodeData) are just ignored.
 * To conclude a successful read of your *.msh file, you should have physical tags to be able to insert elements into a TPZGeoMesh object.
 */
class TPZGmshReader{
    
    /// gmsh file format version (supported versions = {3,4})
    std::string m_format_version;
    
    /// Number of volumes
    int m_n_volumes;
    
    /// Number of surfaces
    int m_n_surfaces;
    
    /// Number of curves
    int m_n_curves;
    
    /// Number of points
    int m_n_points;
    
    /// Number of volumes with physical tag
    int m_n_physical_volumes;
    
    /// Number of surfaces with physical tag
    int m_n_physical_surfaces;
    
    /// Number of curves with physical tag
    int m_n_physical_curves;
    
    /// Number of points with physical tag
    int m_n_physical_points;
    
    /// Geometry dimension
    int m_dimension;
    
    /// Characteristic length to apply a Scale affine transformation
    REAL m_characteristic_lentgh;
    
    //////////// Members related to file format with version 4 ////////////
    
    /// Data structure of both: physical entities dimension and names
    TPZManVector<std::map<int,std::vector<int>>> m_dim_entity_tag_and_physical_tag;
    
    //////////// Members related to file format with version 3 ////////////
    
    /// Structure of both: physical entities dimension and names
    TPZManVector<std::map<int,std::string>,5> m_dim_physical_tag_and_name; /// -> m_dim_physical_tag_and_name
    
    /// Structure of both: names and physical id
    TPZManVector<std::map<std::string,int>,5> m_dim_name_and_physical_tag; /// -> m_dim_name_and_physical_tag
    
    /// Structure of both: dimesion and physical id and "physical tag" @TODO:: Phil please state the need for this
    /// from my point of view this is useless
    TPZManVector<std::map<int,int>,5> m_dim_physical_tag_and_physical_tag; /// m_dim_physical_tag_and_physical_tag
    
    /// Entity index to which the element belongs
    TPZManVector<int64_t> m_entity_index;
    
public:
    
    /// Default constructor
    TPZGmshReader();
    
    /// Default destructor
    ~TPZGmshReader();
    
    /// Copy constructor
    TPZGmshReader(const TPZGmshReader & other);
    
    /// Assignement constructor
    const TPZGmshReader & operator=(const TPZGmshReader & other);
    
    /// Convert Gmsh msh files in a TPZGeoMesh object
    TPZGeoMesh * GeometricGmshMesh(std::string file_name, TPZGeoMesh *gmesh = NULL);

    /// Set the Characteristic length
    void SetCharacteristiclength(REAL length);
    
    /// Set the format version
    void SetFormatVersion(std::string format_version);
    
    /// Print the partition summary after the reading process
    void PrintPartitionSummary(std::ostream & out);
    
    //////////// Members related to file format with version 4 ////////////
    
    /// Convert a Gmsh *.msh file with format 4 to a TPZGeoMesh object
    TPZGeoMesh * GeometricGmshMesh4(std::string file_name, TPZGeoMesh *gmesh = NULL);
    
    void InsertElement(TPZGeoMesh * gmesh, int & physical_identifier, int & el_type, int & el_identifier, std::vector<int> & node_identifiers);
    
    int GetNumberofNodes(int & el_type);
    
    //////////// Members related to file format with version 3 ////////////
    
    /// Convert a Gmsh *.msh file with format 3 to a TPZGeoMesh object
    TPZGeoMesh * GeometricGmshMesh3(std::string file_name, TPZGeoMesh *gmesh = NULL);
    
    /// Insert elements following msh file format */
    bool InsertElement(TPZGeoMesh * gmesh, std::ifstream & line);
    
private:
    
    
};

#endif /* TPZGmshReader_h */
