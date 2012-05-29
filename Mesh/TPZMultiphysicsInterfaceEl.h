/**
 * @file
 * @brief Contains the declaration of multiphysic interface class
 * @author Agnaldo
 * @since 10/26/11.
 */

#ifndef TPZMULTIPHYSICSINTERFACEELH
#define TPZMULTIPHYSICSINTERFACEELH 

#include <iostream>

#include "pzcompel.h"
#include "pzmultiphysicselement.h"


/**
 * @brief Computes the contribution over an interface between two discontinuous elements. \ref CompElement "Computational Element"
 */
class TPZMultiphysicsInterfaceElement : public TPZCompEl {

protected:

	/** @brief Element vector the left of the normal a interface */
	TPZCompElSide 	fLeftElSide;
		
	/** @brief Element vector the right of the normal a interface */
	TPZCompElSide 	fRightElSide;
    
    /** @brief indexes of the connects */
    TPZManVector<int,100> fConnectIndexes;
	
public:
	/** @brief Default constructor */
	TPZMultiphysicsInterfaceElement();	
	/** @brief Constructor */
	TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int &index, TPZCompElSide left, TPZCompElSide right);
    
    /** @brief create a copy of the given element */
    TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy);
    
    /** @brief create a copy of the given element using index mapping */
    TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy, std::map<int,int> & gl2lcConMap,
									std::map<int,int> & gl2lcElMap);
	
	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const 
    {
        return new TPZMultiphysicsInterfaceElement(mesh,*this);
    }
	
	/**
	 * @brief Method for creating a copy of the element in a patch mesh
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the computational elements
     */
	/**
	 * Otherwise of the previous clone function, this method don't
	 * copy entire mesh. Therefore it needs to map the connect index
	 * from the both meshes - original and patch
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
									std::map<int,int> & gl2lcConMap,
									std::map<int,int> & gl2lcElMap) const 
    {
        return new TPZMultiphysicsInterfaceElement(mesh,*this,gl2lcConMap,gl2lcElMap);
    }
	
	/** @brief Default destructor */
	virtual ~TPZMultiphysicsInterfaceElement();
	
	/**
	 * @brief Compute the transform of a paramenter point in the multiphysic interface element to a parameter point in the neighbor super element
	 * @param Neighbor [in] may be this->LeftElementSide() or this->RightElementSide()
	 * @param transf [out] vector of Transforms 
	 */
	void ComputeSideTransform(TPZManVector<TPZCompElSide> &Neighbor, TPZManVector<TPZTransform> &transf);
	
	/**
	 * @brief Maps qsi coordinate at this master element to qsi coordinate at neighbor master element.
	 * @param Neighbor [in] may be this->LeftElementSide() or this->RightElementSide()
	 * @param qsi [in] is the point at this element master
	 * @param NeighIntPoint [out] is the point at neighbor element master. X[qsi] is equal to X[NeighIntPoint]
	 */
	void MapQsi(TPZManVector<TPZCompElSide> &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL> &NeighIntPoint);
    
    /**
     * Add elements to the list of left and right elements
     */
    void SetLeftRightEement(TPZCompElSide &leftel, TPZCompElSide &rightel);
    
    /**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, int index)
    {
        PZError << __FUNCTION__ << " should never be called\n";
        DebugStop();
    }
	
    /** @brief Returns the number of nodes of the element */
	virtual int NConnects() const;
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int ConnectIndex(int i) const;
	

    /** @brief Dimension of the element */
	virtual int Dimension() const
    {
        return Reference()->Dimension();
    }

    /**
     * Compute the stiffness matrix of the interface element
     */
    void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);
    
    /** @brief Initialize the structure of the stiffness matrix */
    void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
    
    /** @brief access function to the left element */
    TPZCompElSide Left() const
    {
        return fLeftElSide;
    }
    
    /** @brief Initialize the material data for the neighbouring element */
    void InitMaterialData(TPZVec<TPZMaterialData> &data, TPZMultiphysicsElement *mfcel);
    
    /** @brief initialize the material data for the geometric data */
    void InitMaterialData(TPZMaterialData &data);
    
    /** @brief Access function to the right element */
    TPZCompElSide Right() const
    {
        return fRightElSide;
    }
    
    virtual int nmeshes()
    {
        if (fLeftElSide) {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(fLeftElSide.Element());
            if(mfcel)
            {
                return mfcel->nmeshes();
            }
        }
        if (fRightElSide) {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element()); 
            if (mfcel) {
                return mfcel->nmeshes();
            }
        }
        return 1;
    }
    
};

#endif