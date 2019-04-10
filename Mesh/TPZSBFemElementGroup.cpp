 //
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzcmesh.h"
#include "tpzintpoints.h"
#include "pzintel.h"
#include "pzelctemp.h"

int TPZSBFemElementGroup::gDefaultPolynomialOrder = 0;



#ifdef USING_MKL
#include <mkl.h>
#elif MACOSX
#include <Accelerate/Accelerate.h>
#endif

#ifdef COMPUTE_CRC
#ifdef USING_BOOST
#include "boost/crc.hpp"
extern TPZVec<boost::crc_32_type::value_type> matglobcrc, eigveccrc, stiffcrc, matEcrc, matEInvcrc;
#endif
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemelementgroup"));
static LoggerPtr loggerMT(Logger::getLogger("pz.mesh.sbfemelementgroupMT"));
static LoggerPtr loggerBF(Logger::getLogger("pz.mesh.sbfemelementgroupBF"));
#endif

#define NONHomogeneous
//#define NONHomogeneous2
#define Prints

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */


TPZSBFemElementGroup::TPZSBFemElementGroup(TPZCompMesh &mesh, int64_t &index) : TPZElementGroup(mesh,index)
{
    fInternalPolynomialOrder = TPZSBFemElementGroup::gDefaultPolynomialOrder;
    
    
//    int nshape = 0;
    int nshape = fInternalPolynomialOrder*(Phi().Rows()) + 1;
    int nvar = 1;
    int64_t newindex = Mesh()->AllocateNewConnect(nshape, nvar, fInternalPolynomialOrder);
    fInternalConnectIndex = newindex;
    
//    int64_t newnodeindex = Mesh()->AllocateNewConnect(nshape, nvar, order);
//    SetConnectIndex(0,newnodeindex);
//    TPZConnect &newnod = Mesh()->ConnectVec()[newnodeindex];
//    int64_t seqnum = newnod.SequenceNumber();
//    Mesh()->Block().Set(seqnum,nvar*nshape);
//    Mesh()->ConnectVec()[ConnectIndex()].IncrementElConnected()
}

void TPZSBFemElementGroup::ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0, TPZElementMatrix &P0, TPZElementMatrix &RF)
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
    InitializeElementMatrix(P0);// LALALALA
    InitializeElementMatrix(RF);// LALALALA
    
    
    
    if(TPZSBFemElementGroup::gDefaultPolynomialOrder){
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
        TPZElementMatrix P0Loc(Mesh(),TPZElementMatrix::EF); //LALALALA
        TPZElementMatrix RFLoc(Mesh(),TPZElementMatrix::EF); //LALALALA
        sbfem->ComputeKMatrices(E0Loc, E1Loc, E2Loc, M0Loc, P0Loc, RFLoc); //LALALAL
        
//        boost::crc_32_type crc;
//        crc.process_bytes(&E0Loc.fMat(0,0),E0Loc.fMat.Rows()*E0Loc.fMat.Rows()*sizeof(STATE));
//        EMatcrc[cel->Index()]= crc.checksum();
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
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
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#ifdef COMPUTE_CRC
        {
            boost::crc_32_type crc;
            int64_t n = E0Loc.fMat.Rows();
            crc.process_bytes(&E0Loc.fMat(0,0), n*n*sizeof(STATE));
            crc.process_bytes(&E1Loc.fMat(0,0), n*n*sizeof(STATE));
            crc.process_bytes(&E2Loc.fMat(0,0), n*n*sizeof(STATE));
            matEcrc[Index()] = crc.checksum();

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
            
            for (int idf = 0; idf<iblsize; idf++) {
                P0.fBlock(ibldest,0,idf,0) += P0Loc.fBlock(ic,0,idf,0); //LALALALA
                RF.fBlock(ibldest,0,idf,0) += RFLoc.fBlock(ic,0,idf,0); //LALALALA
            }
        }
    }
}

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZSBFemElementGroup::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    InitializeElementMatrix(ek, ef);

    // Cálculo de E0, E1, E2, M0 e P0
    if (fComputationMode == EOnlyMass) {
        ek.fMat = fMassMatrix;
        ek.fMat *= fMassDensity;
        ef.fMat.Zero();
        return;
    }
    TPZElementMatrix E0, E1, E2, M0, P0, RF;
    ComputeMatrices(E0, E1, E2, M0, P0, RF);

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
    int dim = Mesh()->Dimension();
    
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
        dgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info); //LU factorization
#endif
#ifdef STATEfloat
        sgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info); //LU factorization
#endif
        if (info != 0) {
            DebugStop();
        }
#ifdef STATEdouble
        dgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info); //Inverse of a LU-factored matrix
#endif
#ifdef STATEfloat
        sgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info); //Inverse of a LU-factored matrix
#endif
        if (info != 0) {
            DebugStop();
        }
    }
    
    // ----- Cálculo de globmat -----
    TPZFMatrix<STATE> globmat(2*n,2*n,0.);
    {
#ifdef STATEdouble
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n); //Matrix-matrix product
#elif defined STATEfloat
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n); //Matrix-matrix product
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
    }
    
    // ----- Solução do problema de autovalor -----
    TPZFMatrix<STATE> globmatkeep(globmat);
    TPZFMatrix< std::complex<double> > eigenVectors;
    TPZManVector<std::complex<double> > eigenvalues;
    
    {
        globmatkeep.SolveEigenProblem(eigenvalues, eigenVectors);
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
    
    // ----- Ordenando autovetores -----
    TPZFNMatrix<200,std::complex<double> > QVectors(n,n,0.);
    fPhi.Resize(n, n); // DÚVIDA: Mudo pro fPhi com todos os autovalores?
    TPZManVector<std::complex<double> > eigvalsel(n,0);
    TPZFMatrix<std::complex<double> > eigvecsel(2*n,n,0.),eigvalmat(1,n,0.);
    
    int count = 0;
    for (int i=0; i<2*n; i++) {
        if (eigenvalues[i].real() < -1.e-6) {
//            double maxvaleigenvec = 0;
            for (int j=0; j<n; j++) {
                QVectors(j,count) = eigenVectors(j+n,i);
                eigvecsel(j,count) = eigenVectors(j,i);
                eigvecsel(j+n,count) = eigenVectors(j+n,i);
                fPhi(j,count) = eigenVectors(j,i);
//                double realvalabs = fabs(fPhi(j,count).real());
//                if (realvalabs > maxvaleigenvec) {
//                    maxvaleigenvec = realvalabs;
//                }
            }
            eigvalsel[count] = eigenvalues[i];
            eigvalmat(0,count) = eigenvalues[i];
//            for (int j=0; j<n; j++) {
//                QVectors(j,count) /= maxvaleigenvec;
//                eigvecsel(j,count) /= maxvaleigenvec;
//                eigvecsel(j+n,count) /= maxvaleigenvec;
//                fPhi(j,count) /= maxvaleigenvec;
//            }
            count++;
        }
    }
    
    TPZFNMatrix<200,std::complex<double> > QVectorsp(n,n,0.); // Adicionei as matrizes relativas aos autovetores positivos
    TPZManVector<std::complex<double> > eigvalselp(n,0);
    TPZFMatrix<std::complex<double> > eigvecselp(2*n,n,0.),eigvalmatp(1,n,0.);
    count=0;
    for (int i=0; i<2*n; i++) {
        if (eigenvalues[i].real() > 1.e-6){
            double maxvaleigenvec = 0;
            for (int j=0; j<n; j++) {
                QVectorsp(j,count) = eigenVectors(j+n,i);
                eigvecselp(j,count) = eigenVectors(j,i);
                eigvecselp(j+n,count) = eigenVectors(j+n,i);
            }
            eigvalselp[count] = eigenvalues[i];
            eigvalmat(0,count) = eigenvalues[i];
            count++;
        }
    }
    
#ifdef PZDEBUG
//    std::cout << "eigval = {" << eigvalsel << "};\n";
//    std::cout << "eigvalpos = {" << eigvalselp << "};\n";
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
                eigvecselp(eq,count) = 1;
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
//    double phinorm = Norm(fPhi);
//    if (loggerMT->isDebugEnabled()) {
//        std::stringstream sout;
//        sout << "Element index " << Index() << " phinorm = " << phinorm;
//        LOGPZ_DEBUG(loggerMT, sout.str())
//    }
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
//    phicopy.Inverse(fPhiInverse, ELU);
//    phicopy.Print("phidec = ",std::cout,EMathematicaInput);
        fPhiInverse.Print("fPhiInverse = ",sout,EMathematicaInput);
//    QVectors.Print("QVectors ", std::cout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZFMatrix<std::complex<double> > ekloc;
    QVectors.Multiply(fPhiInverse, ekloc);
    
    TPZFMatrix<double> ekimag(ekloc.Rows(),ekloc.Cols());
    for (int i=0; i<ekloc.Rows(); i++) {
        for (int j=0; j<ekloc.Cols(); j++) {
            ek.fMat(i,j) = ekloc(i,j).real();
            ekimag(i,j) = ekloc(i,j).imag();
        }
    }
#ifdef COMPUTE_CRC
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
    {
        boost::crc_32_type crc;
        int n = ekloc.Rows();
        crc.process_bytes(&ekloc(0,0), n*n*sizeof(STATE));
        stiffcrc[Index()] = crc.checksum();
    }
    pthread_mutex_unlock(&mutex);
#endif
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
//#ifdef NONHomogeneous
//        sbfem->SetCoefNonHomogeneous(Phi12, A22, A12, fPhi12A22P0, P0.fMat);
//        for (int i=0; i<ef.fMat.Rows(); i++) {
//            for (int j=0; j<ef.fMat.Cols(); j++) {
//                ef.fMat(i,j) = RF.fMat(i,j);
//            }
//        }
//#endif
        
    }
    if (TPZSBFemElementGroup::gDefaultPolynomialOrder) {
        this->CalcStiffBodyLoads(ek, ef);
    } else {
        this->CalcStiff2(ek, ef);
    }
    
//#ifdef NONHomogeneous2
//    this->CollapsedStiffness(E0, E1, E2, P0, ek, ef);
//#endif
//    
//#ifdef NONHomogeneous
//    TPZFMatrix<REAL> partial(n,n,0.), partial2(n,n,0.);
//    TPZFMatrix<REAL> A22Inv(n,n,0.), A12(n,n,0.), A22(n,n,0.);
//    TPZFMatrix<REAL> CoefMatInv(n,n,0.);
//    
//    TPZFMatrix<REAL> Phi11(n,n,0.), Phi21(n,n,0.), Phi12(n,n,0.), Phi22(n,n,0.), Phi11Inv(n,n,0.);
//    {
//        TPZFMatrix<REAL> Phi(2*n, 2*n, 0.);
//        
//        int64_t nrows = eigvecsel.Rows();
//        int64_t ncols = eigvecsel.Cols();
//        for (int64_t i=0; i<nrows; i++) {
//            for (int64_t j=0; j<ncols; j++) {
//                Phi(i,j) = eigvecsel(i,j).real();
//                Phi(i,j+ncols) = eigvecselp(i,j).real();
//            }
//        }
//        for (int64_t i=0; i<n; i++) {
//            for (int64_t j=0; j<n; j++) {
//                Phi11(i,j) = fPhi(i,j).real();
//                Phi11(i,j) = Phi(i,j);
//                Phi21(i,j) = Phi(i+n,j);
//                Phi12(i,j) = Phi(i,j+n);
//                Phi22(i,j) = Phi(i+n,j+n);
//                Phi11Inv(i,j) = fPhiInverse(i,j).real();
//            }
//        }
//        
//        // A22Inv
//        TPZFMatrix<REAL> kstiff = ek.fMat;
//        partial.Zero();
//        kstiff.Multiply(Phi12, partial);
//        A22Inv = Phi22 - partial;
//        nrows = A22Inv.Rows();
//        
//        // A22
//        partial = A22Inv;
//        partial.Resize(n, n-1);
//        partial.PseudoInverse(A22);
//        A22.Resize(n, n);
//        
//        // CoefMatInv
//        int polyorder = 0;
//        REAL coefval = 1 + polyorder + 0.5*dim;
//        for (int64_t i =0; i<n; i++) {
//            if (IsZero(-eigvalselp[i].real() + coefval)) {
//                CoefMatInv(i,i) = 0;
//            }
//            else {
//                CoefMatInv(i,i) = 1/(-eigvalselp[i].real() + coefval);
//            }
//        }
//        
//        // A12
//        partial.Zero();
//        TPZFMatrix<REAL> Phi12Copy(Phi12);
//        Phi12Copy.Multiply(A22, partial);
//        Phi11Inv.Multiply(partial, A12);
//        fA12.Resize(A12.Rows(), A12.Cols());
//        for (int64_t i=0; i<A12.Rows(); i++) {
//            for (int64_t j=0; j<A12.Cols(); j++) {
//                fA12(i,j) = -A12(i,j);
//            }
//        }
//        
//        // fPhi12A22P0
//        fPhi12A22P0.Resize(n, 1);
//        partial.Zero();
//        
//        A22.Multiply(P0.fMat, partial);
//        CoefMatInv.Multiply(partial, partial2);
//        Phi12.Multiply(partial2, fPhi12A22P0);
//    }
//#endif
    
}

void TPZSBFemElementGroup::CalcStiff2(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    
    TPZElementMatrix E0, E1, E2, M0, P0, RF;
    ComputeMatrices(E0, E1, E2, M0, P0, RF);
    
    int n = E0.fMat.Rows();
    
    TPZFNMatrix<100,std::complex<REAL>> K0(n,n,0);
    TPZFNMatrix<100,std::complex<REAL>> rot(n,n,0);
    
    TPZManVector<std::complex<REAL>> eigval(n,0);
    for (int i=0; i<8; i++) {
        eigval[i] = -EigenValues()[i];
    }
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int k=0; k<8; k++) {
                for (int l=0; l<8; l++) {
                    K0(i,j) += (Phi()(k,i) * E0.fMat(k,l) * Phi()(l,j) * (eigval[i]*eigval[j])/(eigval[i]+eigval[j]));
                    K0(i,j) += (Phi()(k,i) * E1.fMat(l,k) * Phi()(l,j) * (eigval[i])/(eigval[i]+eigval[j]));
                    K0(i,j) += (Phi()(k,i) * E1.fMat(k,l) * Phi()(l,j) * (eigval[j])/(eigval[i]+eigval[j]));
                    K0(i,j) += (Phi()(k,i) * E2.fMat(k,l) * Phi()(l,j) / (eigval[i]+eigval[j]));
                }
            }
            rot(i,j) = fPhiInverse(i,j);
        }
    }
    K0(7,7)=0;
    
    TPZFNMatrix<100,std::complex<REAL>> partial(n,n,0);
    K0.Multiply(rot, partial);
    rot.Transpose();
    TPZFNMatrix<100,std::complex<REAL>> K(n,n,0);
    rot.Multiply(partial, K);
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            ek.fMat(i,j) = K(i,j).real();
        }
    }
    
}

void TPZSBFemElementGroup::CalcStiffBodyLoads(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    
    TPZElementMatrix E0, E1, E2, M0, P0, RF;
    ComputeMatrices(E0, E1, E2, M0, P0, RF);
    
    int n = E0.fMat.Rows();
    int norder = this->GetgOrder();
    
    TPZFNMatrix<100,std::complex<REAL>> K0(n + norder*n + 1,n+norder*n + 1,0);
    TPZFNMatrix<100,std::complex<REAL>> rot(n+norder*n + 1, norder*n + 1,0);
    
    TPZManVector<std::complex<double>> eigval(n + norder*n + 1,0);
    for (int i=0; i<n; i++) {
        eigval[i] = -EigenValues()[i];
    }
    if (norder>1) {
        for (int i=0; i<n; i++) {
            eigval[i+n] = 1;
        }
        for (int j=2; j<=norder; j++) {
            for (int i=0; i<n; i++) {
                eigval[i + n*j + 1] = j;
            }
        }
    }
    std::cout << eigval << std::endl;
    
    TPZFNMatrix<200,std::complex<REAL>> Phiu(n,n+n*norder+1);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            Phiu(i,j) = Phi()(i,j);
        }
    }
    if (norder>1) {
        for (int i=0; i<n; i++) {
            Phiu(i,n+i) = 1;
        }
        for (int j=2; j<=norder; j++) {
            for (int i=0; i<n; i++) {
                Phiu(i,n*j+i+1) = 1;
            }
        }
    }
    for (int iorder=0; iorder<norder-1; iorder++) {
        for (int i=0; i<n; i++) {
            Phiu(i,n*(iorder+2)) = Phiu(i,n-1);
        }
    }
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int k=0; k<n; k++) {
                for (int l=0; l<n; l++) {
                    K0(i,j) += (Phiu(k,i).real() * E0.fMat(k,l) * Phiu(l,j) * (eigval[i]*eigval[j])/(eigval[i]+eigval[j]));
                    K0(i,j) += (Phiu(k,i).real() * E1.fMat(l,k) * Phiu(l,j) * (eigval[i])/(eigval[i]+eigval[j]));
                    K0(i,j) += (Phiu(k,i).real() * E1.fMat(k,l) * Phiu(l,j) * (eigval[j])/(eigval[i]+eigval[j]));
                    K0(i,j) += (Phiu(k,i).real() * E2.fMat(k,l) * Phiu(l,j) / (eigval[i]+eigval[j]));
                }
                
                K0(i,j+n) += Phiu(k, i).real() * E0.fMat(k,j) * (eigval[i] * eigval[j + n]) / (eigval[i] + eigval[j + n]);
                K0(i,j+n) += Phiu(k, i).real() * E1.fMat(j,k) * (eigval[i]) / (eigval[i] + eigval[j + n]);
                K0(i,j+n) += Phiu(k, i).real() * E1.fMat(k,j) * (eigval[j+n]) / (eigval[i] + eigval[j + n]);
                K0(i,j+n) += Phiu(k, i).real() * E2.fMat(k,j) / (eigval[i] + eigval[j + n]);
                
                K0(i,j+norder*n+1) += Phiu(k, i).real() * E0.fMat(k,j) * (eigval[i] * eigval[j + norder*n+1]) / (eigval[i] + eigval[j + norder*n+1]);
                K0(i,j+norder*n+1) += Phiu(k, i).real() * E1.fMat(j,k) * (eigval[i]) / (eigval[i] + eigval[j + norder*n+1]);
                K0(i,j+norder*n+1) += Phiu(k, i).real() * E1.fMat(k,j) * (eigval[j + norder*n+1]) / (eigval[i] + eigval[j + norder*n+1]);
                K0(i,j+norder*n+1) += Phiu(k, i).real() * E2.fMat(k,j) / (eigval[i] + eigval[j + norder*n+1]);
                
                K0(i+n,j) += Phiu(k,j).real() * E0.fMat(i,k) * (eigval[i + n] * eigval[j])/(eigval[i + n] + eigval[j]);
                K0(i+n,j) += Phiu(k,j).real() * E1.fMat(k,i) * (eigval[i + n]) / (eigval[i + n] + eigval[j]);
                K0(i+n,j) += Phiu(k,j).real() * E1.fMat(i,k) * (eigval[j])/(eigval[i + n] + eigval[j]);
                K0(i+n,j) += Phiu(k,j).real() * E2.fMat(i,k) / (eigval[i + n] + eigval[j]);
                
                K0(i+norder*n+1,j) += Phiu(k, j).real() * E0.fMat(i,k) * (eigval[i + norder*n+1] * eigval[j])/(eigval[i + norder*n+1] + eigval[j]);
                K0(i+norder*n+1,j) += Phiu(k, j).real() * E1.fMat(k,i) * (eigval[i + norder*n+1])/(eigval[i + norder*n+1] + eigval[j]);
                K0(i+norder*n+1,j) += Phiu(k, j).real() * E1.fMat(i,k) * (eigval[j])/(eigval[i + norder*n+1] + eigval[j]);
                K0(i+norder*n+1,j) += Phiu(k, j).real() * E2.fMat(i,k) / (eigval[i + norder*n+1] + eigval[j]);
                
            }
            
            K0(i+n,j+n) = (E0.fMat(i,j) * eigval[i + n] * eigval[j + n] + E1.fMat(j,i) * eigval[i + n] + E1.fMat(i,j) * eigval[j + n] + E2.fMat(i,j) )/(eigval[i+n] + eigval[j + n]);
            K0(i+norder*n+1,j+norder*n+1) = (E0.fMat(i,j) * eigval[i + norder*n+1] * eigval[j + norder*n+1] + E1.fMat(j,i) * eigval[i + norder*n+1] + E1.fMat(i,j) * eigval[j + norder*n+1] + E2.fMat(i,j))/(eigval[i + norder*n+1] + eigval[j + norder*n+1]);
            K0(i+n,j+norder*n+1) = (E0.fMat(i,j) * eigval[i + n] * eigval[j + norder*n+1] + E1.fMat(j,i) * eigval[i + n] + E1.fMat(i,j) * eigval[j + norder*n+1] + E2.fMat(i,j)) / (eigval[i + n] + eigval[j + norder*n+1]);
            K0(i+norder*n+1,j+n) = (E0.fMat(i,j) * eigval[i + norder*n+1] * eigval[j + n] + E1.fMat(j,i) * eigval[i + norder*n+1] + E1.fMat(i,j) * eigval[j + n] + E2.fMat(i,j))/(eigval[i + norder*n+1] + eigval[j + n]);
            K0(i+norder*n+1,j+n) = (E0.fMat(i,j) * eigval[i + norder*n+1] * eigval[j + n] + E1.fMat(j,i) * eigval[i + norder*n+1] + E1.fMat(i,j) * eigval[j + n] + E2.fMat(i,j))/(eigval[i + norder*n+1] + eigval[j + n]);
            
            rot(i,j) = fPhiInverse(i,j);
        }
    }
    K0(n-1,n-1)=0;
    
    for (int i=0; i<n; i++) {
        rot(i+n, i+n+1) = 4;
        rot(i+norder*n+1, i+n+1) = -4;
    }
    
    for (int i=0; i<n; i++) {
        rot(i+n,n) = -Phi()(i,7).real();
    }
    rot(norder*n,n)=1;
    
    TPZFNMatrix<100,std::complex<REAL>> partial(n+norder*n+1,norder*n+1,0);
    K0.Multiply(rot, partial);
    rot.Transpose();
    TPZFNMatrix<100,std::complex<REAL>> K(norder*n+1,norder*n+1,0);
    rot.Multiply(partial, K);
    
    ek.fMat.Resize(norder*n+1,norder*n+1);
    for (int i=0; i<norder*n+1; i++) {
        for (int j=0; j<norder*n+1; j++) {
            ek.fMat(i,j) = K(i,j).real();
        }
    }
    
    // force vector
    int64_t nel = fElGroup.size();
    TPZFNMatrix<200,std::complex<double>> f(n+norder*n+1,1,0);
    
    for (int i=0; i<n+norder*n+1; i++) {
        for (int64_t j = 0; j<nel; j++) {
            TPZCompEl *cel = fElGroup[j];
            TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
            if (!sbfem) {
                DebugStop();
            }
#endif
            sbfem->LocalBodyForces(f, Phiu, eigval, i, j);
        }
    }
//    for (int i=0; i<EigenValues().size()+1; i++) {
//        fEigenvalues[i] = eigval[i];
//    }
    
    ef.fMat.Resize(n*norder+1,1);
    ef.fMat.Zero();
    TPZFNMatrix<200,std::complex<double>> efcomplex;
    rot.Multiply(f, efcomplex);
    for (int i=0; i<n*norder+1; i++) {
        ef.fMat(i,0) = efcomplex(i,0).real();
    }
    fek = ek.fMat;
    fef = ef.fMat;
    
    std::cout << ek.fSourceIndex << std::endl;
    
#ifdef LOG4CXX
    if (loggerBF->isDebugEnabled()) {
        std::stringstream sout;
        E0.fMat.Print("E0 = ", sout, EMathematicaInput);
        E1.fMat.Print("E1 = ", sout, EMathematicaInput);
        E2.fMat.Print("E2 = ", sout, EMathematicaInput);
        Phiu.Print("Phiu = ", sout, EMathematicaInput);
        rot.Transpose();
        rot.Print("rot = ", sout, EMathematicaInput);
        ek.fMat.Print("ek = ", sout, EMathematicaInput);
        ef.fMat.Print("ef = ", sout, EMathematicaInput);
        f.Print("f = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerBF, sout.str())
    }
#endif
    
    fEigenvalues = eigval;
    fPhi.Zero();
    Phiu.Multiply(rot, fPhi);
    
    for (int i=0; i<n+norder*n+1; i++) {
        for (int64_t j = 0; j<nel; j++) {
            TPZCompEl *cel = fElGroup[j];
            TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
            if (!sbfem) {
                DebugStop();
            }
#endif
            sbfem->SetPhiEigVal(fPhi, fEigenvalues);
        }
    }
    
}

void TPZSBFemElementGroup::LoadSolution()
{
    int nc = NConnects();
    int ndof = fPhi.Rows();
    if(TPZSBFemElementGroup::gDefaultPolynomialOrder){
        ndof=0;
        for (int64_t ic=0; ic<nc; ic++) {
            ndof += Mesh()->ConnectVec()[fConnectIndexes[ic]].NShape()* Mesh()->ConnectVec()[fConnectIndexes[ic]].NState();
        }
    }
    
    TPZFNMatrix<200,std::complex<double> > sol(ndof, fMesh->Solution().Cols(),0.);
    fCoef.Resize(ndof, fPhi.Cols());
    TPZFNMatrix<200,std::complex<double> > unh_local(ndof, fMesh->Solution().Cols(),0.);
    
    int count = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = Connect(ic);
        int nshape = c.NShape();
        int nstate = c.NState();
        int blsize = nshape*nstate;
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = fMesh->Block().Position(seqnum);
        for (int seq=0; seq < blsize; seq++) {
            for (int c=0; c<sol.Cols(); c++)
            {
                sol(count+seq,c) = fMesh->Solution()(pos+seq,c);
//#ifdef NONHomogeneous
//                unh_local(count+seq,c) = -fPhi12A22P0(count+seq,0);
//#endif
            }
        }
        count += blsize;
    }
    if(TPZSBFemElementGroup::gDefaultPolynomialOrder){
//        fCoef = Mesh()
    } else{
        fPhiInverse.Multiply(sol-unh_local, fCoef);
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
    int64_t ncon = NConnects();
    int64_t nshape = 0;
    
    std::vector<int64_t> connects(ncon);
    for (int64_t i=0; i<ncon; i++) {
        connects[i] = fConnectIndexes[i];
        nshape += Mesh()->ConnectVec()[fConnectIndexes[i]].NShape();
    }
    
    std::vector<int64_t>::iterator it = std::find(connects.begin(), connects.end(), fInternalConnectIndex);
    if (it != connects.end()) {
        return;
    }
    else{
        fConnectIndexes.resize(ncon+1);
        fConnectIndexes[ncon] = fInternalConnectIndex;
        
        TPZConnect &c = Mesh()->ConnectVec()[fInternalConnectIndex];
        c.SetNShape(nshape);
        int64_t seq = c.SequenceNumber();
        Mesh()->Block().Set(seq, nshape);
    }
}



//http://www.netlib.org/lapack/lug/node50.html
//https://software.intel.com/en-us/node/521079
