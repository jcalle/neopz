/**
 * @file pzcondensedcompel.cpp
 * @brief Contains the implementations of the TPZElementGroup methods.
 */

#include "pzelementgroup.h"
#include "pzlog.h"
#include "pzstepsolver.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzelementgroup"));
#endif

TPZElementGroup::TPZElementGroup() : TPZCompEl(), fElGroup(), fConnectIndexes()
{
}


TPZElementGroup::~TPZElementGroup()
{
    if (fElGroup.size()) {
        DebugStop();
    }
}

/** @brief create a copy of the condensed computational element in the other mesh */
TPZElementGroup::TPZElementGroup(TPZCompMesh &mesh, const TPZElementGroup &copy) : TPZCompEl(mesh, copy)
{
    TPZStack<TPZCompEl *> newel;
    int nel = copy.fElGroup.size();
    for (int el=0; el<nel; el++) {
        newel.Push(copy.fElGroup[el]->Clone(mesh) );
    }
    for (int el=0; el<nel; el++) {
        AddElement(newel[el]);
    }
}

/** @brief add an element to the element group
 */
void TPZElementGroup::AddElement(TPZCompEl *cel)
{
    fElGroup.Push(cel);
    std::set<int> connects;
    int nc = fConnectIndexes.size();
    for (int ic=0; ic<nc; ic++) {
        connects.insert(fConnectIndexes[ic]);
    }
    nc = cel->NConnects();
    for (int ic=0; ic<nc; ic++) {
        connects.insert(cel->ConnectIndex(ic));
    }
    nc = connects.size();
    if (nc != fConnectIndexes.size()) {
        fConnectIndexes.Resize(nc, 0);
        std::set<int>::iterator it = connects.begin();
        for (int ic = 0; it != connects.end(); it++,ic++) {
            fConnectIndexes[ic] = *it;
        }
    }
    int elindex = cel->Index();
    Mesh()->ElementVec()[elindex] = 0;
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Hiding element index " << elindex << " from the mesh data structure";
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
}

/** @brief unwrap the condensed element from the computational element and delete the condensed element */
void TPZElementGroup::Unwrap()
{
    int nel = fElGroup.size();
    for (int el=0; el<nel; el++) {
        int elindex = fElGroup[el]->Index();
        Mesh()->ElementVec()[elindex] = fElGroup[el];
    }
    fElGroup.Resize(0);
    fConnectIndexes.Resize(0);
    delete this;
}

/**
 * @brief Set the index i to node inode
 * @param inode node to set index
 * @param index index to be seted
 */
void TPZElementGroup::SetConnectIndex(int inode, int index)
{
    LOGPZ_ERROR(logger,"SetConnectIndex should never be called")
    DebugStop();
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
TPZCompEl *TPZElementGroup::ClonePatchEl(TPZCompMesh &mesh,
                                std::map<int,int> & gl2lcConMap,
                                std::map<int,int> & gl2lcElMap) const
{
    int index;
    TPZElementGroup *result = new TPZElementGroup(mesh,index);
    int nel = fElGroup.size();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el]->ClonePatchEl(mesh,gl2lcConMap,gl2lcElMap);
        result->AddElement(cel);
    }
    return result;
}

void TPZElementGroup::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) const {
    InitializeElementMatrix(ek);
    int rows = ek.fMat.Rows();
    ek.fMat.Redim(rows, rows);
    ek.fType = TPZElementMatrix::EK;
    InitializeElementMatrix(ef);
}//void

void TPZElementGroup::InitializeElementMatrix(TPZElementMatrix &ef) const {
	const int ncon = this->NConnects();
	int numeq = 0;
	ef.fBlock.SetNBlocks(ncon);
    for (int ic=0; ic<ncon; ic++) {
        TPZConnect &c = Connect(ic);
        int blsize = c.NShape()*c.NState();
        numeq += blsize;
        ef.fBlock.Set(ic, blsize);
    }
    const int numloadcases = 1;
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
	ef.fMat.Redim(numeq,numloadcases);
	ef.fNumStateVars = 1;
	ef.fConnect.Resize(ncon);
	for(int i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
	}
}//void


/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZElementGroup::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    std::map<int,int> locindex;
    int ncon = fConnectIndexes.size();
    for (int ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
    InitializeElementMatrix(ek, ef);
    int nel = fElGroup.size();
    TPZElementMatrix ekloc,efloc;
    for (int el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        cel->CalcStiff(ekloc, efloc);
        int nelcon = ekloc.NConnects();
        for (int ic=0; ic<nelcon; ic++) {
            int iblsize = ekloc.fBlock.Size(ic);
            int icindex = ekloc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int idf = 0; idf<iblsize; idf++) {
                ef.fBlock(ibldest,0,idf,0) += efloc.fBlock(ic,0,idf,0);
            }
            for (int jc = 0; jc<nelcon; jc++) {
                int jblsize = ekloc.fBlock.Size(jc);
                int jcindex = ekloc.fConnect[jc];
                int jbldest = locindex[jcindex];
                for (int idf = 0; idf<iblsize; idf++) {
                    for (int jdf=0; jdf<jblsize; jdf++) {
                        ek.fBlock(ibldest,jbldest,idf,jdf) += ekloc.fBlock(ic,jc,idf,jdf);
                    }
                }
            }
        }        
    }
}


/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
void TPZElementGroup::CalcResidual(TPZElementMatrix &ef)
{
    std::map<int,int> locindex;
    int ncon = fConnectIndexes.size();
    for (int ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
    InitializeElementMatrix(ef);
    int nel = fElGroup.size();
    TPZElementMatrix ekloc,efloc;
    for (int el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        cel->CalcStiff(ekloc, efloc);
        int nelcon = ekloc.NConnects();
        for (int ic=0; ic<nelcon; ic++) {
            int iblsize = ekloc.fBlock.Size(ic);
            int icindex = ekloc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int idf = 0; idf<iblsize; idf++) {
                ef.fBlock(ibldest,idf,0,0) += efloc.fBlock(ic,idf,0,0);
            }
        }        
    }
}

