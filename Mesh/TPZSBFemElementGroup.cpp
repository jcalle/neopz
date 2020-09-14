//
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#include "pzelementgroup.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzcmesh.h"
#include "tpzintpoints.h"
#include "pzintel.h"
#include "pzelctemp.h"
#include "TPZMaterial.h"

#include <numeric>

int TPZSBFemElementGroup::gDefaultPolynomialOrder = 0;

int TPZSBFemElementGroup::gPolynomialShapeFunctions = 0;

#ifdef USING_MKL
#include <mkl.h>
#elif MACOSX
#include <Accelerate/Accelerate.h>
#endif

#ifdef COMPUTE_CRC
#ifdef USING_BOOST
#include "boost/crc.hpp"
extern TPZVec<boost::crc_32_type::value_type> matglobcrc, eigveccrc, stiffcrc, matEcrc, matEInvcrc, matPhicrc;
#endif
#endif

#ifdef USING_BLAZE
//#define BLAZE_DEFAULT_ALIGMENT_FLAG = blaze::unaligned
//#define BLAZE_DEFAULT_PADDING_FLAG = blaze::unpadded
#define BLAZE_USE_VECTORIZATION 0
#define BLAZE_USE_SHARED_MEMORY_PARALLELIZATION 0
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/config/Thresholds.h>
#include <blaze/Math.h>
using blaze::columnMajor;
using blaze::DynamicMatrix;
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemelementgroup"));
static LoggerPtr loggercoefmatrices(Logger::getLogger("pz.mesh.sbfemcoefmatrices"));
static LoggerPtr loggerMT(Logger::getLogger("pz.mesh.sbfemelementgroupMT"));
static LoggerPtr loggerBF(Logger::getLogger("pz.mesh.sbfemelementgroupBF"));
static LoggerPtr loggerbubble(Logger::getLogger("pz.mesh.sbfembubbleparam"));
#endif

TPZSBFemElementGroup::TPZSBFemElementGroup(TPZCompMesh &mesh, int64_t &index) : TPZElementGroup(mesh,index)
{
    fInternalPolynomialOrder = TPZSBFemElementGroup::gDefaultPolynomialOrder;
    if (fInternalPolynomialOrder != 0) {
        int nshape = fInternalPolynomialOrder*(Phi().Rows());
        int nvar = 1;
        int64_t newindex = Mesh()->AllocateNewConnect(nshape, nvar, fInternalPolynomialOrder);
        fInternalConnectIndex = newindex;
    }
}

/**
 * @brief Computes the SBFEM coefficient matrices and the mass matrix
 * @param E0 E0 coefficient matrix
 * @param E1 E1 coefficient matrix
 * @param E2 E2 coefficient matrix
 * @param M0 Mass matrix
 */
void TPZSBFemElementGroup::ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0)
{
    std::map<int64_t,int64_t> locindex;
    int64_t ncon = fConnectIndexes.size();
    for (int64_t ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
    TPZElementMatrix ef(Mesh(),TPZElementMatrix::EF);
    int64_t nel = fElGroup.size();
    InitializeElementMatrix(E0, ef);
    InitializeElementMatrix(E1, ef);
    InitializeElementMatrix(E2, ef);
    InitializeElementMatrix(M0, ef);

    if(TPZSBFemElementGroup::gDefaultPolynomialOrder != 0){
        int ndof=0;
        for (int64_t ic=0; ic<ncon-1; ic++) {
            ndof += Mesh()->ConnectVec()[fConnectIndexes[ic]].NShape()* Mesh()->ConnectVec()[fConnectIndexes[ic]].NState();
        }
        E0.fMat.Resize(ndof, ndof);
        E1.fMat.Resize(ndof, ndof);
        E2.fMat.Resize(ndof, ndof);
    }

    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        TPZElementMatrix E0Loc(Mesh(),TPZElementMatrix::EK);
        TPZElementMatrix E1Loc(Mesh(),TPZElementMatrix::EK);
        TPZElementMatrix E2Loc(Mesh(),TPZElementMatrix::EK);
        TPZElementMatrix M0Loc(Mesh(),TPZElementMatrix::EK);
        sbfem->ComputeKMatrices(E0Loc, E1Loc, E2Loc,M0Loc);
        
#ifdef LOG4CXX
        if (loggercoefmatrices->isDebugEnabled()) {
            TPZGeoEl *gel = cel->Reference();
            
            int matid = 0;
            if(gel) matid = gel->MaterialId();
            std::stringstream sout;
            if (gel) {
                sout << "Material id " << matid <<std::endl;
            }
            else
            {
                sout << "No associated geometry\n";
            }
            sout << "Connect indexes ";
            for (int i=0; i<cel->NConnects(); i++) {
                sout << cel->ConnectIndex(i) << " ";
            }
            sout << std::endl;
            sout << "Local indexes ";
            for (int i=0; i<cel->NConnects(); i++) {
                sout << locindex[cel->ConnectIndex(i)] << " ";
            }
            sout << std::endl;
            E0Loc.fMat.Print("Matriz elementar E0",sout);
            E1Loc.fMat.Print("Matriz elementar E1",sout);
            E2Loc.fMat.Print("Matriz elementar E2",sout);
            M0Loc.fMat.Print("Matriz elementar M0",sout);
            LOGPZ_DEBUG(loggercoefmatrices, sout.str())
        }
#endif
        int nelcon = E0Loc.NConnects();
        for (int ic=0; ic<nelcon; ic++) {
            int iblsize = E0Loc.fBlock.Size(ic);
            int icindex = E0Loc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int jc = 0; jc<nelcon; jc++) {
                int jblsize = E0Loc.fBlock.Size(jc);
                int jcindex = E0Loc.fConnect[jc];
                int jbldest = locindex[jcindex];
                for (int idf = 0; idf<iblsize; idf++) {
                    for (int jdf=0; jdf<jblsize; jdf++) {
                        E0.fBlock(ibldest,jbldest,idf,jdf) += E0Loc.fBlock(ic,jc,idf,jdf);
                        E1.fBlock(ibldest,jbldest,idf,jdf) += E1Loc.fBlock(ic,jc,idf,jdf);
                        E2.fBlock(ibldest,jbldest,idf,jdf) += E2Loc.fBlock(ic,jc,idf,jdf);
                        M0.fBlock(ibldest,jbldest,idf,jdf) += M0Loc.fBlock(ic,jc,idf,jdf);
                    }
                }
            }
        }
    }
}

/**
 * @brief Compute the eigenvalues and eigenvectors using blaze-lib
 */
void TPZSBFemElementGroup::ComputeEigenvaluesBlaze()
{
#ifdef USING_BLAZE
    
    TPZElementMatrix E0, E1, E2, M0;
    ComputeMatrices(E0, E1, E2, M0);

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "BLAZE VERSION\n";
        E0.fMat.Print("E0 = ",sout, EMathematicaInput);
        E1.fMat.Print("E1 = ",sout, EMathematicaInput);
        E2.fMat.Print("E2 = ",sout, EMathematicaInput);
        M0.fMat.Print("M0 = ",sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    int n = E0.fMat.Rows();
    
    int dim = Mesh()->Dimension();
    
    blaze::DynamicMatrix<STATE,blaze::columnMajor> E0blaze(n,n), E1blaze(n,n), E2blaze(n,n), E0Invblaze(n,n);

    memcpy(&E0blaze.data()[0], E0.fMat.Adress(), E0blaze.spacing()*n*sizeof(STATE));
    memcpy(&E1blaze.data()[0], E1.fMat.Adress(), E0blaze.spacing()*n*sizeof(STATE));
    memcpy(&E2blaze.data()[0], E2.fMat.Adress(), E0blaze.spacing()*n*sizeof(STATE));

#ifdef COMPUTE_CRC
    {
        boost::crc_32_type crc;
        if(E0blaze.spacing() != n) DebugStop();
        crc.process_bytes(&E0blaze(0,0), E0blaze.spacing()*n*sizeof(STATE));
        crc.process_bytes(&E1blaze(0,0), E1blaze.spacing()*n*sizeof(STATE));
        crc.process_bytes(&E2blaze(0,0), E2blaze.spacing()*n*sizeof(STATE));
        matEcrc[Index()] = crc.checksum();
    }
#endif
    
    E0Invblaze = inv( E0blaze );  // Compute the inverse of E0

#ifdef COMPUTE_CRC
    {
        boost::crc_32_type crc;
        if(E0Invblaze.spacing() != n) DebugStop();
        crc.process_bytes(&E0Invblaze(0,0), E0Invblaze.spacing()*n*sizeof(STATE));
        matEInvcrc[Index()] = crc.checksum();
    }
#endif

    blaze::DynamicMatrix<STATE,blaze::columnMajor> globmatblaze(2*n,2*n);
    blaze::DynamicMatrix<STATE,blaze::columnMajor> E0InvE1Tblaze = E0Invblaze*trans(E1blaze);
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
        	globmatblaze(i, j) = E0InvE1Tblaze(i,j);
            globmatblaze(i, j+n) = -E0Invblaze(i,j);
        }
    }

    blaze::DynamicMatrix<STATE,blaze::columnMajor> E1E0InvE1T = E1blaze*E0InvE1Tblaze;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmatblaze(i+n,j) = E1E0InvE1T(i,j)-E2blaze(i,j);
        }
    }

    blaze::DynamicMatrix<STATE,blaze::columnMajor> E1E0Invblaze = E1blaze*E0Invblaze;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++){
            globmatblaze(i+n,j+n) = -E1E0Invblaze(i,j);
        }
    }
    for(int i=0; i<n; i++){
    	globmatblaze(i,i) -= (dim-2)*0.5;
        globmatblaze(i+n,i+n) += (dim-2)*0.5;
    }
    
    blaze::DynamicVector<blaze::complex<double>,blaze::columnVector> eigvalblaze( 2*n );       // The vector for the real eigenvalues
	blaze::DynamicMatrix<blaze::complex<double>,blaze::columnMajor> eigvecblaze( 2*n, 2*n );
    
    // The matrix for the left eigenvectors
    static pthread_mutex_t mutex_serial = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex_serial);
	eigen(globmatblaze, eigvalblaze, eigvecblaze);
    pthread_mutex_unlock(&mutex_serial);

#ifdef COMPUTE_CRC2
    pthread_mutex_lock(&mutex);
    extern int gnumthreads;
    std::stringstream sout;
    sout << "eigval" << gnumthreads << ".nb";
    static int count_loc = 0;
    std::ofstream file;
    if (count_loc == 0) {
        file.open(sout.str());
    }
    else
    {
        file.open(sout.str(),std::ios::app);
    }
    std::stringstream eigv;
    eigv << "EigVec" << Index() << " = ";
    if(count_loc < 1)
    {
        eigvalblaze.Print(eigv.str().c_str(),file,EMathematicaInput);
    }
    count_loc++;
    pthread_mutex_unlock(&mutex);
#endif
    
    fQVectors.Resize(n, n);
    fPhi.Resize(n, n);
    TPZManVector<std::complex<double> > eigvalsel(n,0);
    TPZFMatrix<std::complex<double> > eigvecsel(2*n,n,0.),eigvalmat(1,n,0.);
    int count = 0;
    for (int i=0; i<2*n; i++) {
        if (eigvalblaze[i].real() < -1.e-6) {
            double maxvaleigenvec = 0;
            for (int j=0; j<n; j++) {
                fQVectors(j,count) = eigvecblaze(j+n,i);
                eigvecsel(j,count) = eigvecblaze(j,i);
                eigvecsel(j+n,count) = eigvecblaze(j+n,i);
                fPhi(j,count) = eigvecblaze(j,i);
                double realvalabs = fabs(fPhi(j,count).real());
                if (realvalabs > maxvaleigenvec) {
                    maxvaleigenvec = realvalabs;
                }
            }
            eigvalsel[count] = eigvalblaze[i];
            eigvalmat(0,count) = eigvalblaze[i];
            for (int j=0; j<n; j++) {
                fQVectors(j,count) /= maxvaleigenvec;
                eigvecsel(j,count) /= maxvaleigenvec;
                eigvecsel(j+n,count) /= maxvaleigenvec;
                fPhi(j,count) /= maxvaleigenvec;
            }
            count++;
        }
    }
    
#ifdef PZDEBUG
   std::cout << "eigval = {" << eigvalsel << "};\n";
#endif

    if (dim == 2)
    {
        int nstate = Connect(0).NState();
        if (nstate != 2 && nstate != 1) {
            DebugStop();
        }
        if(count != n-nstate) {
            DebugStop();
        }
        int ncon = fConnectIndexes.size();
        int eq=0;
        std::set<int64_t> cornercon;
        BuildCornerConnectList(cornercon);
        for (int ic=0; ic<ncon; ic++) {
            int64_t conindex = ConnectIndex(ic);
            if (cornercon.find(conindex) != cornercon.end())
            {
                fPhi(eq,count) = 1;
                eigvecsel(eq,count) = 1;
                if (nstate == 2)
                {
                    fPhi(eq+1,count+1) = 1;
                    eigvecsel(eq+1,count+1) = 1;
                }
            }
            eq += Connect(ic).NShape()*Connect(ic).NState();
        }
    }
    if(dim==3 && count != n)
    {
        std::cout << __PRETTY_FUNCTION__ << __LINE__ << " count = " << count << " n = " << n << std::endl;
        for(int i=0; i< 2*n; i++) std::cout << eigvalblaze[i] << std::endl;
        DebugStop();
    }
    fEigenvalues = eigvalsel;
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "eigenvalues BLAZE " << eigvalsel << std::endl;
        fPhi.Print("Phivec =",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef COMPUTE_CRC
    {
        boost::crc_32_type crc;
        crc.process_bytes(&fPhi(0,0), fPhi.Rows()*fPhi.Cols()*sizeof(STATE));
        matPhicrc[Index()] = crc.checksum();
    }
#endif
        
    fPhiInverse.Redim(n, n);
    blaze::DynamicMatrix<std::complex<double>,blaze::columnMajor> phiblaze(n,n), PhiInverseblaze(n,n);
    memcpy(&phiblaze.data()[0], fPhi.Adress(),n*n*sizeof(std::complex<double>));
    PhiInverseblaze = inv( phiblaze );  // Compute the inverse of A
    if (PhiInverseblaze.spacing() != n) DebugStop();
    memcpy(fPhiInverse.Adress(), &PhiInverseblaze.data()[0], PhiInverseblaze.spacing()*n*sizeof(std::complex<double>));


#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        fPhiInverse.Print("fPhiInverse = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif

    if(1)
    {
        std::ofstream out("EigenProblem.nb");
        TPZFMatrix<STATE> globmatkeep(2*n,2*n,0);
        memcpy(globmatkeep.Adress(), &globmatblaze.data()[0], 2*n*2*n*sizeof(STATE));
        globmatkeep.Print("matrix = ",out,EMathematicaInput);
        eigvecsel.Print("eigvec =",out,EMathematicaInput);
        eigvalmat.Print("lambda =",out,EMathematicaInput);
        fPhi.Print("phi = ",out,EMathematicaInput);
        fPhiInverse.Print("phiinv = ",out,EMathematicaInput);
        fQVectors.Print("qvec = ",out,EMathematicaInput);
    }
    
    ComputeMassMatrix(M0);

#ifdef COMPUTE_CRC
    static pthread_mutex_t mutex =PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex);
    {
        boost::crc_32_type crc;
        crc.process_bytes(&eigvecblaze(0,0), n*n*sizeof(std::complex<double>));
        eigveccrc[Index()] = crc.checksum();
    }
    {
        boost::crc_32_type crc;
        int n = ekloc.Rows();
        crc.process_bytes(&ekloc(0,0), n*n*sizeof(STATE));
        stiffcrc[Index()] = crc.checksum();
    }
    {
        boost::crc_32_type crc;
        crc.process_bytes(&globmatblaze(0,0), globmatblaze.spacing()*n*sizeof(STATE));
        matglobcrc[Index()] = crc.checksum();
    }
    pthread_mutex_unlock(&mutex);
#endif
#endif
}

/**
 * @brief Compute the eigenvalues and eigenvectors using Intel MKL
 */
void TPZSBFemElementGroup::ComputeEigenvaluesMKL()
{

    if (fElGroup.size() == 0) return;

#ifdef USING_BLAZE
    CalcStiffBlaze(ek,ef);
    return;
#endif
    InitializeElementMatrix(ek, ef);

    if (fComputationMode == EOnlyMass) {
        ek.fMat = fMassMatrix;
        ek.fMat *= fMassDensity;
        ef.fMat.Zero();
        return;
    }

    TPZElementMatrix E0,E1,E2, M0;
    ComputeMatrices(E0, E1, E2, M0);

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        E0.fMat.Print("E0 = ",sout, EMathematicaInput);
        E1.fMat.Print("E1 = ",sout, EMathematicaInput);
        E2.fMat.Print("E2 = ",sout, EMathematicaInput);
        M0.fMat.Print("M0 = ",sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    int n = E0.fMat.Rows();
    
    int dim = this->Mesh()->Dimension();
    
    TPZFMatrix<STATE> E0Inv(E0.fMat);
    if(0)
    {
        try
        {
            TPZFMatrix<STATE> E0copy(E0.fMat);
            TPZVec<int> pivot;
            E0copy.Decompose_LU(pivot);
            E0Inv.Identity();
            E0copy.Substitution(&E0Inv, pivot);
        }
        catch(...)
        {
            exit(-1);
        }
    }
    else
    {
        TPZVec<int> pivot(E0Inv.Rows(),0);
        int nwork = 4*n*n + 2*n;
        TPZVec<STATE> work(2*nwork,0.);
        int info=0;
#ifdef STATEdouble
        dgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
#endif
#ifdef STATEfloat
        sgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
#endif
        if (info != 0) {
            DebugStop();
        }
#ifdef STATEdouble
        dgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
#endif
#ifdef STATEfloat
        sgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
#endif
        if (info != 0) {
            DebugStop();
        }
    }
    
    TPZFMatrix<STATE> globmat(2*n,2*n,0.);
    
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i,j+n) = -E0Inv(i,j);
        }
    }
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1.fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1.fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i+n,j) -= E2.fMat(i,j);
        }
    }

#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1.fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1.fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif

    for (int i=0; i<n; i++) {
        globmat(i,i) -= (dim-2)*0.5;
        globmat(i+n,i+n) += (dim-2)*0.5;
    }
    
    static pthread_mutex_t mutex_serial = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex_serial);
    
    TPZFMatrix<STATE> globmatkeep(globmat);
    TPZFMatrix< std::complex<double> > eigenVectors;
    TPZManVector<std::complex<double> > eigenvalues;
    globmatkeep.SolveEigenProblem(eigenvalues, eigenVectors);
    
    pthread_mutex_unlock(&mutex_serial);
 
    {
#ifdef COMPUTE_CRC
        static pthread_mutex_t mutex =PTHREAD_MUTEX_INITIALIZER;
        pthread_mutex_lock(&mutex);
        extern int gnumthreads;
        std::stringstream sout;
        sout << "eigval" << gnumthreads << ".nb";
        static int count = 0;
        std::ofstream file;
        if (count == 0) {
            file.open(sout.str());
        }
        else
        {
            file.open(sout.str(),std::ios::app);
        }
        std::stringstream eigv;
        eigv << "EigVec" << Index() << " = ";
        if(count < 1)
        {
            eigenVectors.Print(eigv.str().c_str(),file,EMathematicaInput);
        }
        count++;
        pthread_mutex_unlock(&mutex);
#endif
    }

    if(0)
    {
        TPZManVector<STATE> eigvalreal(2*n,0.);
        TPZFMatrix<STATE> eigvecreal(2*n,2*n,0.);
        for (int i=0; i<2*n; i++) {
            eigvalreal[i] = eigenvalues[i].real();
            for (int j=0; j<2*n; j++) {
                eigvecreal(i,j) = eigenVectors(i,j).real();
            }
        }
    }
        
    fQVectors.Resize(n, n);
    fPhi.Resize(n, n);
    TPZManVector<std::complex<double> > eigvalsel(n,0);
    TPZFMatrix<std::complex<double> > eigvecsel(2*n,n,0.),eigvalmat(1,n,0.);
    int count = 0;
    for (int i=0; i<2*n; i++) {
        if (eigenvalues[i].real() < -1.e-6) {
            double maxvaleigenvec = 0;
            for (int j=0; j<n; j++) {
                fQVectors(j,count) = eigenVectors(j+n,i);
                eigvecsel(j,count) = eigenVectors(j,i);
                eigvecsel(j+n,count) = eigenVectors(j+n,i);
                fPhi(j,count) = eigenVectors(j,i);
                double realvalabs = fabs(fPhi(j,count).real());
                if (realvalabs > maxvaleigenvec) {
                    maxvaleigenvec = realvalabs;
                }
            }
            eigvalsel[count] = eigenvalues[i];
            eigvalmat(0,count) = eigenvalues[i];
            for (int j=0; j<n; j++) {
                fQVectors(j,count) /= maxvaleigenvec;
                eigvecsel(j,count) /= maxvaleigenvec;
                eigvecsel(j+n,count) /= maxvaleigenvec;
                fPhi(j,count) /= maxvaleigenvec;

            }
            count++;
        }
    }
    
#ifdef PZDEBUG
    std::cout << "eigval = {" << eigvalsel << "};\n";
#endif

    if (dim == 2)
    {
        int nstate = Connect(0).NState();
        if (nstate != 2 && nstate != 1) {
            DebugStop();
        }
        if(count != n-nstate) {
            DebugStop();
        }
        int ncon = fConnectIndexes.size();
        int eq=0;
        std::set<int64_t> cornercon;
        BuildCornerConnectList(cornercon);
        for (int ic=0; ic<ncon; ic++) {
            int64_t conindex = ConnectIndex(ic);
            if (cornercon.find(conindex) != cornercon.end())
            {
                fPhi(eq,count) = 1;
                eigvecsel(eq,count) = 1;
                if (nstate == 2)
                {
                    fPhi(eq+1,count+1) = 1;
                    eigvecsel(eq+1,count+1) = 1;
                }
            }
            eq += Connect(ic).NShape()*Connect(ic).NState();
        }
    }
    if(dim==3 && count != n)
    {
        DebugStop();
    }

    fEigenvalues = eigvalsel;

#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "eigenvalues " << eigvalsel << std::endl;
        fPhi.Print("Phivec =",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
   double phinorm = Norm(fPhi);
   if (loggerMT->isDebugEnabled()) {
       std::stringstream sout;
       sout << "Element index " << Index() << " phinorm = " << phinorm;
       LOGPZ_DEBUG(loggerMT, sout.str())
   }
#endif
        
    TPZFMatrix<std::complex<double> > phicopy(fPhi);
    fPhiInverse.Redim(n, n);
    fPhiInverse.Identity();
    
    try
    {
        TPZVec<int> pivot;
        phicopy.Decompose_LU(pivot);
        phicopy.Substitution(&fPhiInverse, pivot);
    }
    catch(...)
    {
        exit(-1);
    }

#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        phicopy.Inverse(fPhiInverse, ELU);
        phicopy.Print("phidec = ", sout,EMathematicaInput);
        fPhiInverse.Print("fPhiInverse = ",sout,EMathematicaInput);
        fQVectors.Print("QVectors ", sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

#ifdef COMPUTE_CRC
    static pthread_mutex_t mutex =PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex);
    {
        boost::crc_32_type crc;
        int64_t n = E0.fMat.Rows();
        crc.process_bytes(&E0.fMat(0,0), n*n*sizeof(STATE));
        crc.process_bytes(&E1.fMat(0,0), n*n*sizeof(STATE));
        crc.process_bytes(&E2.fMat(0,0), n*n*sizeof(STATE));
        matEcrc[Index()] = crc.checksum();
        
    }
    {
        boost::crc_32_type crc;
        int64_t n = E0Inv.Rows();
        crc.process_bytes(&E0Inv(0,0), n*n*sizeof(STATE));
        matEInvcrc[Index()] = crc.checksum();
        
    }
    {
        boost::crc_32_type crc;
        crc.process_bytes(&globmat(0,0), 4*n*n*sizeof(STATE));
        matglobcrc[Index()] = crc.checksum();
    }
    {
        boost::crc_32_type crc;
        crc.process_bytes(&eigenVectors(0,0), n*n*sizeof(std::complex<double>));
        eigveccrc[Index()] = crc.checksum();
    }
    pthread_mutex_unlock(&mutex);
#endif
    if(1)
    {
        std::ofstream out("EigenProblem.nb");
        eigvecsel.Print("eigvec =",out,EMathematicaInput);
        eigvalmat.Print("lambda =",out,EMathematicaInput);
        fPhi.Print("phi = ",out,EMathematicaInput);
        fPhiInverse.Print("phiinv = ",out,EMathematicaInput);
        fQVectors.Print("qvec = ",out,EMathematicaInput);
    }

}

/**
 * @brief Compute the eigenvalues and the bubble parameters
 */
void TPZSBFemElementGroup::ComputeEigenvalues()
{
#ifdef USING_BLAZE
    ComputeEigenvaluesBlaze();
#else
    ComputeEigenvaluesMKL();
#endif

    if (fInternalPolynomialOrder > 0){
        ComputeBubbleParameters();
    }

    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->SetPhiEigVal(fPhi, fEigenvalues);
        sbfem->SetCoefNonHomogeneous(fPhiBubble, fEigenvaluesBubble, fPhiInverse, fMatBubble);
    }

}

void TPZSBFemElementGroup::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    InitializeElementMatrix(ek, ef);

    TPZElementMatrix E0,E1,E2, M0;
    ComputeMatrices(E0, E1, E2, M0);

    if (fComputationMode == EOnlyMass) {
        ek.fMat = fMassMatrix;
        ek.fMat *= fMassDensity;
        ef.fMat.Zero();
        return;
    }

    TPZFMatrix<std::complex<double> > ekloc;
    fQVectors.Multiply(fPhiInverse, ekloc);
    
    TPZFMatrix<double> ekimag(ekloc.Rows(),ekloc.Cols());
    for (int i=0; i<ekloc.Rows(); i++) {
        for (int j=0; j<ekloc.Cols(); j++) {
            ek.fMat(i,j) = ekloc(i,j).real();
            ekimag(i,j) = ekloc(i,j).imag();
        }
    }
    
    ComputeMassMatrix(M0);
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        ek.fMat.Print("Stiff = ",sout,EMathematicaInput);
        fMassMatrix.Print("Mass = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    if(fComputationMode == EMass)
    {
        int nr = ek.fMat.Rows();
        for (int r=0; r<nr; r++) {
            for (int c=0; c<nr; c++) {
                ek.fMat(r,c) += fMassMatrix(r,c)/fDelt;
            }
        }
    }

    if (fInternalPolynomialOrder > 0)
    {
        // Computing the stiffness matrix related to the bubbles
        int64_t nbubbles = fEigenvaluesBubble.size();
        int n = fPhi.Rows();

        TPZFNMatrix<100, std::complex<REAL>> *k01 = reinterpret_cast<TPZFNMatrix<100, std::complex<REAL>>*>(&(E0.fMat));
        TPZFNMatrix<100, std::complex<REAL>> *k02 = reinterpret_cast<TPZFNMatrix<100, std::complex<REAL>>*>(&(E1.fMat));
        TPZFNMatrix<100, std::complex<REAL>> *k03 = reinterpret_cast<TPZFNMatrix<100, std::complex<REAL>>*>(&(E2.fMat));

        TPZFNMatrix<100,std::complex<REAL>> temp1(n,nbubbles,0), temp2(n,nbubbles,0), temp3(n,nbubbles,0);

        k01->Multiply(fPhiBubble, temp1); 
        k02->Multiply(fPhiBubble, temp2);
        k03->Multiply(fPhiBubble, temp3);

        bool transpose = 1;
        fPhiBubble.Multiply(temp1, *k01, transpose); // k01 = fPhi^T * E0 * fPhi
        fPhiBubble.Multiply(temp2, *k02, transpose); // k02 = fPhi^T * E1 * fPhi
        fPhiBubble.Multiply(temp3, *k03, transpose); // k02 = fPhi^T * E2 * fPhi
        
        TPZFNMatrix<200,std::complex<REAL>> K0(nbubbles,nbubbles,0);
        // K0 = ( k01* (-eigval[i])*(-eigval[j]) + k02^T * (-eigval[i]) + k02 * (-eigval[j]) + k03 / (-eigval[i]-eigval[j])
        for (int i=0; i<nbubbles; i++) {
            for (int j=0; j<nbubbles; j++) {
                if((-fEigenvaluesBubble[i]-fEigenvaluesBubble[j]).real() == 0) {
                    K0(i,j) += 0;
                } else {
                    K0(i,j) = (k01->GetVal(i,j) * (-fEigenvaluesBubble[i]*-fEigenvaluesBubble[j]) + k02->GetVal(j,i) * -fEigenvaluesBubble[i] +
                     k02->GetVal(i,j) * -fEigenvaluesBubble[j] + k03->GetVal(i,j))/(-fEigenvaluesBubble[i]-fEigenvaluesBubble[j]);
                }
            }
        }
    
        K0.Multiply(fMatBubble, temp1);
        
        int ndofbubbles = fMatBubble.Cols();
        TPZFNMatrix<100,std::complex<REAL>> K(ndofbubbles,ndofbubbles,0);
        fMatBubble.Multiply(temp1, K, transpose);  // K = fMatBubble^T K fMatBubble;
        
        ek.fMat.Resize(n+ndofbubbles,n+ndofbubbles);
        for (int i=0; i<ndofbubbles; i++) {
            for (int j=0; j<ndofbubbles; j++) {
                ek.fMat(i+n,j+n) = K(i,j).real();
            }
        }
        
#ifdef LOG4CXX
        if(loggerbubble->isDebugEnabled())
        {
            std::stringstream sout;
            K.Print("KBubble = ",sout,EMathematicaInput);
            LOGPZ_DEBUG(loggerbubble, sout.str())
        }
#endif

        // Computing the force vector
        int icon = this->ConnectIndex(NConnects()-1);

        TPZFNMatrix<200,std::complex<REAL>> f(n,1,0);
        TPZFNMatrix<200,std::complex<REAL>> fbubble(nbubbles,1,0);

        int64_t nel = fElGroup.size();
        
        for (int64_t j = 0; j<nel; j++) {
            TPZCompEl *cel = fElGroup[j];
            TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
            if (!sbfem) {
                DebugStop();
            }
#endif
            sbfem->LocalBodyForces(f, fbubble, fEigenvalues, fEigenvaluesBubble, icon);
        }
        ef.fMat.Zero();

        TPZFNMatrix<200,std::complex<double>> efcomplex;
        TPZFNMatrix<200,std::complex<double>> efcomplexbubbles;

        fPhiInverse.Multiply(f, efcomplex, transpose);
        fMatBubble.Multiply(fbubble, efcomplexbubbles, transpose);
        
        ef.fMat.Resize(n+ndofbubbles,1);
        for (int i=0; i<n; i++) {
            ef.fMat(i,0) = -efcomplex(i,0).real();
        }
        for (int i=0; i<ndofbubbles; i++) {
            ef.fMat(i+n,0) = -efcomplexbubbles(i,0).real();
        }
        
#ifdef LOG4CXX
        if (loggerBF->isDebugEnabled()) {
            std::stringstream sout;

            K0.Print("K0 = ", sout, EMathematicaInput);
            f.Print("f = ", sout, EMathematicaInput);
            fbubble.Print("fbubble = ", sout, EMathematicaInput);
            
            ek.fMat.Print("ek = ", sout, EMathematicaInput);
            ef.fMat.Print("ef = ", sout, EMathematicaInput);
            LOGPZ_DEBUG(loggerBF, sout.str())
        }
#endif
    }
}

void TPZSBFemElementGroup::LoadSolution()
{
    
    int nc = NConnects();
    int ndof = fPhiInverse.Cols()+fMatBubble.Cols();
    
    TPZFNMatrix<200,std::complex<double> > uh_local(ndof, fMesh->Solution().Cols(),0.);
    fCoef.Resize(ndof,fMesh->Solution().Cols());
    
    int count = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = Connect(ic);
        int nshape = c.NShape();
        int nstate = c.NState();
        int blsize = nshape*nstate;
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = fMesh->Block().Position(seqnum);
        for (int seq=0; seq < blsize; seq++) {
            for (int c=0; c<uh_local.Cols(); c++)
            {
                uh_local(count+seq,c) = fMesh->Solution()(pos+seq,c);
            }
        }
        count += blsize;
    }
    if(fInternalPolynomialOrder > 0){
        TPZFMatrix<std::complex<double> > coef0, coefbubble;
        TPZFMatrix<std::complex<double> > uh_local0(uh_local), uh_localbubble(fMatBubble.Rows(),0);
        uh_local0.Resize(fPhiInverse.Rows(),1);
        for (int64_t i = 0; i < fMatBubble.Rows(); i++)
        {
            uh_localbubble(i,0) = uh_local(i+fPhiInverse.Rows());
        }
        fPhiInverse.Multiply(uh_local0, coef0);
        fMatBubble.Multiply(uh_localbubble, coefbubble);
        fCoef.Resize(ndof,1);
        for (int64_t i = 0; i < fPhiInverse.Cols(); i++)
        {
            fCoef(i,0) = coef0(i);
        }
        for (int64_t i = 0; i < fMatBubble.Cols(); i++)
        {
            fCoef(i+fPhiInverse.Cols(),0) = coefbubble(i);
        }
    } else{
        fPhiInverse.Multiply(uh_local, fCoef);
    }

    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->LoadCoef(fCoef);
    }
    #ifdef LOG4CXX
        if (loggerBF->isDebugEnabled()) {
            std::stringstream sout;
            fCoef.Print("fCoef = ", sout, EMathematicaInput);
            LOGPZ_DEBUG(loggerBF, sout.str())
        }
#endif

}

/// Compute the mass matrix based on the value of M0 and the eigenvectors
void TPZSBFemElementGroup::ComputeMassMatrix(TPZElementMatrix &M0)
{
    //    M0.fMat.Print("Mass = ",std::cout,EMathematicaInput);
    TPZFMatrix<std::complex<double> > temp;
    REAL alpha = 1.;
    REAL beta = 0.;
    int transpose = 1;
    int nrow = fEigenvalues.size();
    TPZFMatrix<std::complex<double> > M0Loc(nrow,nrow),MassLoc(nrow,nrow);
    for (int i=0; i<nrow; i++) {
        for(int j=0; j<nrow; j++)
        {
            M0Loc(i, j) = M0.fMat(i,j);
        }
    }
    fPhi.MultAdd(M0Loc, M0Loc, temp, alpha, beta, transpose);
    //    temp.Print("Temp = ",std::cout,EMathematicaInput);
    temp.Multiply(fPhi,MassLoc);
    for (int i=0; i<nrow; i++) {
        for (int j=0; j<nrow; j++) {
            MassLoc(i,j) /= (2.-fEigenvalues[i]-fEigenvalues[j]);
        }
    }
    fPhiInverse.MultAdd(MassLoc, M0Loc, temp,alpha, beta, transpose);
    temp.Multiply(fPhiInverse,MassLoc);
    
    TPZFMatrix<double> MassLocImag(nrow,nrow);
    fMassMatrix.Redim(nrow, nrow);
    for (int i=0; i<nrow; i++) {
        for(int j=0; j<nrow; j++)
        {
            fMassMatrix(i, j) = MassLoc(i,j).real();
            MassLocImag(i,j) = MassLoc(i,j).imag();
        }
    }
    
//    fMassMatrix.Print("Mass Matrix",std::cout,EMathematicaInput);
#ifdef PZDEBUG
//    std::cout << "Norm of imaginary part " << Norm(MassLocImag) << std::endl;
#endif
}

void TPZSBFemElementGroup::LoadEigenVector(int64_t eig)
{
    fCoef.Zero();
    fCoef(eig,0) = 1.;
    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->LoadCoef(fCoef);
    }
    
}

/** @brief add an element to the element group
 */
void TPZSBFemElementGroup::AddElement(TPZCompEl *cel)
{
    std::set<int64_t> connects;
    int nc = fConnectIndexes.size();
    for (int ic=0; ic<nc; ic++) {
        connects.insert(fConnectIndexes[ic]);
    }
    TPZSBFemVolume *celvol = dynamic_cast<TPZSBFemVolume *>(cel);
    TPZCompEl *celskeleton = Mesh()->Element(celvol->SkeletonIndex());
    nc = celskeleton->NConnects();
    for (int ic=0; ic<nc; ic++) {
        connects.insert(celskeleton->ConnectIndex(ic));
    }
    nc = connects.size();
    if (nc != fConnectIndexes.size()) {
        fConnectIndexes.Resize(nc, 0);
        std::set<int64_t>::iterator it = connects.begin();
        for (int ic = 0; it != connects.end(); it++,ic++) {
            fConnectIndexes[ic] = *it;
        }
    }
    TPZElementGroup::AddElement(cel);
}

void TPZSBFemElementGroup::InitializeInternalConnect()
{
    if(fElGroup.size() == 0) return;
    int64_t ncon = NConnects();
    int nshapeboundary = 0;
    int64_t nshapeinternal = 0;

    std::vector<int64_t> connects(ncon);
    for (int64_t i=0; i<ncon; i++) {
        connects[i] = fConnectIndexes[i];
        nshapeboundary += Mesh()->ConnectVec()[fConnectIndexes[i]].NShape()*Mesh()->ConnectVec()[fConnectIndexes[i]].NState();
    }
    if(!(fElGroup[0]->Material()) ){
        std::cout << "TPZSBFemElementGroup::InitializeInternalConnect - No material associated\n";
        DebugStop();
    }
    nshapeinternal += (fInternalPolynomialOrder-1)*nshapeboundary; // if p=1, only the hat function is added
    int nstate = fElGroup[0]->Material()->NStateVariables();
    nshapeinternal += nstate; // hat function
    
    // Verifies if the internal connect is allocated into the connect list
    std::vector<int64_t>::iterator it = std::find(connects.begin(), connects.end(), fInternalConnectIndex);
    if (it != connects.end()) {
        return;
    }
    else{
        fConnectIndexes.resize(ncon+1);
        fConnectIndexes[ncon] = fInternalConnectIndex;
        
        TPZConnect &c = Mesh()->ConnectVec()[fInternalConnectIndex];
        c.SetNShape(nshapeinternal);
        int64_t seq = c.SequenceNumber();
        Mesh()->Block().Set(seq, nshapeinternal);
    }
}

void TPZSBFemElementGroup::ComputeBubbleParameters()
{
    int cont = 0;
    int nstate = fElGroup[0]->Material()->NStateVariables();
    int n = fPhi.Rows();

    // Finding the rational eigenvalues
    TPZManVector<int64_t> ind(0);
    for (int i=0; i<n; i++) {
        REAL resto = fEigenvalues[i].real() - nearbyint(fEigenvalues[i].real());
        if( !IsZero(resto)){
            ind.resize(ind.size()+1);
            ind[cont] = i;
            cont++;
        }
    }
    fEigenvaluesBubble.resize(cont*2 + n*fInternalPolynomialOrder + nstate);
    for (int i = 0; i < cont; ++i)
    {
        fEigenvaluesBubble[2*i] = fEigenvalues[ind[i]].real();
        fEigenvaluesBubble[2*i+1] = -fInternalPolynomialOrder;
    }
    int64_t neigval = fEigenvaluesBubble.size();
    
    // Creating new connects for the bubble functions
    TPZConnect &c = Mesh()->ConnectVec()[fInternalConnectIndex];
    int64_t seq = c.SequenceNumber();
    int internalpoly = Mesh()->Block().Size(seq);
    int neq_adc = internalpoly + cont;
    c.SetNShape(neq_adc);
    Mesh()->Block().Set(seq, neq_adc);
    Mesh()->ExpandSolution();
    
    // Updating the eigenvalue vector for polynomial bubbles
    for (int i=0; i<n; i++) {
        fEigenvaluesBubble[i + 2*cont] = -1;
    }
    for (int j=0; j<nstate; j++) {
        fEigenvaluesBubble[n + 2*cont + j] = 0;
    }
    for (int j=2; j<=fInternalPolynomialOrder; j++) {
        for (int i=0; i<n; i++) {
            fEigenvaluesBubble[i + n*(j-1) + nstate + 2*cont] = -j;
        }
    }

    // Updating the eigenvectors matrix
    int64_t nphis = fEigenvaluesBubble.size();
    fPhiBubble.Resize(n, nphis);
    fPhiBubble.Zero();
    for (int i=0; i<n; i++) {
        for (int j=0; j<cont; j++) {
            fPhiBubble(i,j*2) = fPhi(i,ind[j]).real();
            fPhiBubble(i,j*2+1) = -fPhi(i,ind[j]).real();
        }
    }
    for (int i=0; i<n; i++) {
        fPhiBubble(i,i+2*cont) = 1;
        for (int j=0; j<nstate; j++) {
            fPhiBubble(i,n+2*cont+j) = fPhi(i,n-nstate+j);
        }
    }
    for (int j=1; j<fInternalPolynomialOrder; j++) {
        for (int i=0; i<n; i++) {
            fPhiBubble(i,n*j+i + 2*cont + nstate) = 1;
        }
    }
    
    // Updating the fPhiInverse
    fMatBubble.Resize(fInternalPolynomialOrder*n + nstate + cont*2, (fInternalPolynomialOrder-1)*n + nstate + cont);
    
    // rational functions
    for (int64_t i=0; i<cont; i++) {
        fMatBubble(2*i, i) = 1;
        fMatBubble(2*i+1, i) = 1;
    }
    // hat functions
    for (int64_t i=0; i<n; i++) {
        for (int j=0; j<nstate; j++) {
            fMatBubble(i + 2*cont, j + cont) = -1.;
            if(IsZero(fPhi(i, n - nstate + j).real())){
                fMatBubble(2*cont + i, j + cont) = 0.;
            }
        }
    }
    //polynomial functions
    for (int j=0; j<nstate; j++) {
        fMatBubble(n + j + 2*cont,j + cont) = 1;
    }
    if(fInternalPolynomialOrder > 1){
        for (int i=0; i<n; i++) {
            fMatBubble(i+ 2*cont, i+1 + cont+ nstate-1) = 1;
            fMatBubble(n+1 + i + 2*cont+ nstate-1, i+1 + cont+ nstate-1) = -1;
        }
    }
    for (int j=2; j<fInternalPolynomialOrder; j++) {
        for (int i=0; i<n; i++) {
            fMatBubble(i+n*(j-1) + 2*cont + nstate, i+n*(j-1) + cont+ nstate) = 1;
            fMatBubble(i+n*j + 2*cont+ nstate, i+n*(j-1) + cont+ nstate) = -1;
        }
    }

#ifdef LOG4CXX
    if (loggerbubble->isDebugEnabled()) {
        std::stringstream sout;
        sout << "eigvalbubbles = {" << fEigenvaluesBubble << "};\n";
        fPhiBubble.Print("fPhiBubble = ", sout, EMathematicaInput);
        fMatBubble.Print("fMatBubble = ", sout, EMathematicaInput);
        // fMat.Print("fMat = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerbubble, sout.str())
    }
#endif

}

void TPZSBFemElementGroup::OverwritePhis(TPZFMatrix<std::complex<double>> &Phiu, TPZFNMatrix<200,std::complex<double>> &matphiinv, TPZManVector<std::complex<double> > &eigval)
{
    int64_t n = Phiu.Rows();
    
    for (int i=0; i<n; i++) {
        matphiinv(i+n,n) = -Phiu(i,n-1).real();
    }
    matphiinv(2*n,n) = 1;
    
    for (int64_t i=0; i<n; i++) {
        for (int64_t j=0; j<n; j++) {
            Phiu(i,j) = 0.;
            matphiinv(i,j) = 0.;
        }
        Phiu(i,i) = 1.;
        matphiinv(i,i) = 1.;
        eigval[i] = -1;
    }
}

//http://www.netlib.org/lapack/lug/node50.html
//https://software.intel.com/en-us/node/521079
