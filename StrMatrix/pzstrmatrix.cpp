// -*- c++ -*-

//$Id: pzstrmatrix.cpp,v 1.28 2010-02-18 20:32:55 phil Exp $

/* Generated by Together */

#include "pzstrmatrix.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"

#include "pzgnode.h"
#include "TPZTimer.h"
#include "TPZMTAssemble.h"

#include "pzcheckconsistency.h"

using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.tpzstructmatrix"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

TPZStructMatrix::TPZStructMatrix(TPZCompMesh *mesh) : fMinEq(-1), fMaxEq(-1) {
    fMesh = mesh;
}

TPZStructMatrix::~TPZStructMatrix() {}

TPZMatrix *TPZStructMatrix::Create() {
  cout << "TPZStructMatrix::Create should never be called\n";
  return 0;
}

TPZStructMatrix *TPZStructMatrix::Clone() {
  cout << "TPZStructMatrix::Clone should never be called\n";
  return 0;
}

void TPZStructMatrix::Assemble(TPZMatrix & stiffness, TPZFMatrix & rhs){
  if (fMesh){
    TPZStructMatrix::Assemble(stiffness, rhs, *fMesh,fMinEq,fMaxEq, &fMaterialIds);
  }
}

void TPZStructMatrix::pthread_Assemble(TPZMatrix & stiffness, TPZFMatrix & rhs){
	if (fMesh){
		TPZStructMatrix::pthread_Assemble(stiffness, rhs, *fMesh,fMinEq,fMaxEq, &fMaterialIds);
	}
}

void TPZStructMatrix::Assemble(TPZFMatrix & rhs){
	if (fMesh){
		TPZStructMatrix::Assemble(rhs, *fMesh,fMinEq,fMaxEq, &fMaterialIds);
	}
	else 
	{
			LOGPZ_ERROR(logger,"Assemble called without mesh")
	}

}

void TPZStructMatrix::pthread_Assemble(TPZFMatrix & rhs)
{
	if (fMesh){
		TPZStructMatrix::pthread_Assemble(rhs, *fMesh,fMinEq,fMaxEq, &fMaterialIds);
	}
	else 
	{
		LOGPZ_ERROR(logger,"pthread_Assemble called without mesh")
	}

}

/** STATIC METHOD **/
#ifdef WIN32
#include "threadExecuteProdIndex.h"

void TPZStructMatrix::Assemble(TPZMatrix & stiffness, TPZFMatrix & rhs, TPZCompMesh &mesh,
															int mineq, int maxeq, std::set<int> *MaterialIds,
															TExecuteProdIndex* myThread){
#else
void TPZStructMatrix::Assemble(TPZMatrix & stiffness, TPZFMatrix & rhs, TPZCompMesh &mesh,
															int mineq, int maxeq, std::set<int> *MaterialIds){
#endif
  int iel;
  int nelem = mesh.NElements();
  TPZElementMatrix ek(&mesh, TPZElementMatrix::EK),ef(&mesh, TPZElementMatrix::EF);
	bool globalresult = true;
	bool writereadresult = true;

#ifndef WIN32
	TPZTimer calcstiff("Computing the stiffness matrices");
	TPZTimer assemble("Assembling the stiffness matrices");
#endif
  TPZAdmChunkVector<TPZCompEl *> &elementvec = mesh.ElementVec();

	int count = 0;
  for(iel=0; iel < nelem; iel++) {
    TPZCompEl *el = elementvec[iel];
    if(!el) continue;

    if(MaterialIds){
      TPZAutoPointer<TPZMaterial> mat = el->Material();
      if (!mat) continue;
      int matid = mat->Id();
      if (MaterialIds->find(matid) == MaterialIds->end()) continue;
    }///if

	  count++;
	  if(!(count%20))
	  {
		  std::cout << '*';
		  std::cout.flush();
	  }
	  if(!(count%500))
	  {
		  std::cout << "\n";
	  }
#ifndef WIN32
		calcstiff.start();
#endif

#ifdef WIN32
	if(myThread) if(myThread->AmIKilled()){
		return;
	}
#endif

		el->CalcStiff(ek,ef);
	  
#ifdef CHECKCONSISTENCY
	  //extern TPZCheckConsistency stiffconsist("ElementStiff");
	  stiffconsist.SetOverWrite(true);
	  bool result;
	  result = stiffconsist.CheckObject(ek.fMat);
	  if(!result)
	  {
		  globalresult = false;
		  std::stringstream sout;
		  sout << "element " << iel << " computed differently";
		  LOGPZ_ERROR(loggerCheck,sout.str())
	  }
	  
	  /*
	  static TPZCheckConsistency checkstiff("ElementStiff");
	  checkstiff.SetReadMode();
	  checkstiff.SetOverWrite(true);
	  result = checkstiff.CheckObject(ek.fMat);
	  if(!result)
	  {
		  writereadresult = false;
		  std::stringstream sout;
		  sout << "element " << iel << " computed differently";
		  LOGPZ_ERROR(loggerCheck,sout.str())
	  }
	  */
#endif
	  
#ifndef WIN32
		calcstiff.stop();
		assemble.start();
#endif

    if(!el->HasDependency()) {
      ek.ComputeDestinationIndices();
      if(mineq != -1 || maxeq != -1)
      {
        FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,mineq,maxeq);
      }
      stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
      rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
#ifdef LOG4CXX
		if(loggerel->isDebugEnabled())
		{
			std::stringstream sout;
			ek.fMat.Print("Element stiffness matrix",sout);
			ef.fMat.Print("Element right hand side", sout);
			LOGPZ_DEBUG(loggerel,sout.str())
		}
#endif
    } else {
      // the element has dependent nodes
      ek.ApplyConstraints();
      ef.ApplyConstraints();
      ek.ComputeDestinationIndices();
      if(mineq != -1 || maxeq != -1)
      {
        FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,mineq,maxeq);
      }
      stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
      rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
		}

#ifndef WIN32
		assemble.stop();
#endif		
  }//fim for iel
	std::cout << std::endl;
	{
		std::stringstream sout;
		sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
		LOGPZ_DEBUG(loggerCheck,sout.str())
	}
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << "Number of equations " << stiffness.Rows() << std::endl;
    sout << calcstiff.processName() << " " << calcstiff << std::endl;
    sout << assemble.processName() << " " << assemble;
/*    stiffness.Print("Matriz de Rigidez: ",sout);
    rhs.Print("Vetor de Carga: ",sout);*/
    LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif

}
	
void TPZStructMatrix::pthread_Assemble(TPZMatrix & mat, TPZFMatrix & rhs,
						  TPZCompMesh &mesh,
						  int mineq,int maxeq,
						  std::set<int> *MaterialIds)
{
	ThreadData threaddata(mesh,mat,rhs,mineq,maxeq,MaterialIds);
	
	int numthreads = 5;
	std::vector<pthread_t> allthreads(numthreads-1);
	int itr;
	for(itr=0; itr<numthreads-1; itr++)
	{
		pthread_create(&allthreads[itr], NULL,ThreadData::ThreadWork, &threaddata);
	}
	
	ThreadData::ThreadAssembly(&threaddata);

	for(itr=0; itr<numthreads-1; itr++)
	{
		pthread_join(allthreads[itr],NULL);
	}
	
}
	
void TPZStructMatrix::pthread_Assemble(TPZFMatrix & rhs, TPZCompMesh &mesh, int mineq, int maxeq, std::set<int> *MaterialIds)
{
	Assemble(rhs, mesh, mineq, maxeq, MaterialIds);
}

/** STATIC METHOD **/
// ofstream effile("serialEF.txt");
// #define MTCALCSTIFF
void TPZStructMatrix::Assemble(TPZFMatrix & rhs, TPZCompMesh &mesh, int mineq, int maxeq, std::set<int> *MaterialIds){

#ifdef MTCALCSTIFF
  TPZMTAssemble::AssembleMT(rhs, mesh, mineq, maxeq, MaterialIds);
  return;
#endif

  int iel;
  int nelem = mesh.NElements();
//   TPZElementMatrix ef(&mesh, TPZElementMatrix::EF);

#ifndef WIN32
	TPZTimer calcresidual("Computing the residual vector");
	TPZTimer assemble("Assembling the residual vector");
#endif

  TPZAdmChunkVector<TPZCompEl *> &elementvec = mesh.ElementVec();

  for(iel=0; iel < nelem; iel++) {
    TPZCompEl *el = elementvec[iel];
    if(!el) continue;

    if(MaterialIds){
      TPZAutoPointer<TPZMaterial> mat = el->Material();
      if (!mat) continue;
      int matid = mat->Id();
      if (MaterialIds->find(matid) == MaterialIds->end()) continue;
    }///if
    
    TPZElementMatrix ef(&mesh, TPZElementMatrix::EF);

#ifndef WIN32
		calcresidual.start();
#endif

		el->CalcResidual(ef);

#ifndef WIN32
		calcresidual.stop();
#endif
    
//     effile << "\n************** " << el->Reference()->Type() << " **************\n";
//     el->Print(effile);
//     ef.Print(effile);
//     effile.flush();

#ifndef WIN32
		assemble.start();
#endif

    if(!el->HasDependency()) {
      ef.ComputeDestinationIndices();
      if(mineq != -1 & maxeq != -1)
      {
        FilterEquations(ef.fSourceIndex,ef.fDestinationIndex,mineq,maxeq);
      }
      rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
    } else {
      // the element has dependent nodes
      ef.ApplyConstraints();
      ef.ComputeDestinationIndices();
      if(mineq != -1 & maxeq != -1)
      {
        FilterEquations(ef.fSourceIndex,ef.fDestinationIndex,mineq,maxeq);
      }
      rhs.AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
		}

#ifndef WIN32
		assemble.stop();
#endif

  }//fim for iel
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << calcresidual.processName() << " " << calcresidual << std::endl;
    sout << assemble.processName() << " " << assemble;
    LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif

}

  /// filter out the equations which are out of the range
void TPZStructMatrix::FilterEquations(TPZVec<int> &origindex, TPZVec<int> &destindex, int mineq, int upeq)
{
  if(mineq == -1 || upeq == -1) return;
  int count = 0;
  int nel = origindex.NElements();
  int i;
  for(i=0; i<nel; i++)
  {
    if(destindex[i] >= mineq && destindex[i] < upeq)
    {
      origindex[count] = origindex[i];
      destindex[count++] = destindex[i]-mineq;
    }
  }
  origindex.Resize(count);
  destindex.Resize(count);
  
}

TPZMatrix * TPZStructMatrix::CreateAssemble(TPZFMatrix &rhs)
{
    TPZMatrix *stiff = Create();
    int neq = stiff->Rows();
    rhs.Redim(neq,1);
    Assemble(*stiff,rhs);
	
#ifdef LOG4CXX
	if(loggerel->isDebugEnabled())
	{
		std::stringstream sout;
		stiff->Print("Stiffness matrix",sout);
		rhs.Print("Right hand side", sout);
		LOGPZ_DEBUG(loggerel,sout.str())
	}
#endif
    return stiff;
	
}
/*
struct ThreadData
	{
		// Initialize the mutex semaphores and others
		ThreadData();
		// Destroy the mutex semaphores and others
		~ThreadData();
		// current structmatrix object
		TPZStructMatrix *fReference;
		// mutexes (to choose which element is next)
		pthread_mutex_t fAccessElement;
		// semaphore (to wake up assembly thread)
		sem_t fAssembly;
		// list of computed element matrices (autopointers?)
		std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > 
		// current element
		int fNextElement;
		// look for an element index which needs to be computed
		int FindElement();
		
	};
*/

	TPZStructMatrix::ThreadData::ThreadData(TPZCompMesh &mesh, TPZMatrix &mat, TPZFMatrix &rhs, int mineq, int maxeq, std::set<int> *MaterialIds) : fNextElement(0), fMesh(&mesh),
	fGlobMatrix(&mat), fGlobRhs(&rhs), fMinEq(mineq), fMaxEq(maxeq)
{
	if(MaterialIds) 
	{
		fMaterialIds = *MaterialIds;
	}
	pthread_mutex_init(&fAccessElement,NULL);
/*	sem_t *sem_open( ... );
    int sem_close(sem_t *sem);
    int sem_unlink(const char *name);
*/	
	std::stringstream sout;
	static int counter = 0;
	sout << "AssemblySem" << counter++;
	fAssembly = sem_open(sout.str().c_str(), O_CREAT,777,1);
	if(fAssembly == SEM_FAILED)
	{
		std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
	}
}

TPZStructMatrix::ThreadData::~ThreadData()
	{
		pthread_mutex_destroy(&fAccessElement);
		sem_close(fAssembly);
	}

void *TPZStructMatrix::ThreadData::ThreadWork(void *datavoid)
{
	ThreadData *data = (ThreadData *) datavoid;
	// compute the next element (this method is threadsafe)
	int iel = data->NextElement();
	int mineq = data->fMinEq;
	int maxeq = data->fMaxEq;
	TPZCompMesh *cmesh = data->fMesh;
	int nel = cmesh->NElements();
	while(iel < nel)
	{
		TPZAutoPointer<TPZElementMatrix> ek = new TPZElementMatrix(cmesh,TPZElementMatrix::EK);
		TPZAutoPointer<TPZElementMatrix> ef = new TPZElementMatrix(cmesh,TPZElementMatrix::EF);

		TPZCompEl *el = cmesh->ElementVec()[iel];
		TPZElementMatrix *ekp = ek.operator->();
		TPZElementMatrix *efp = ef.operator->();
		TPZElementMatrix &ekr = *ekp;
		TPZElementMatrix &efr = *efp;
		el->CalcStiff(ekr,efr);
		
		if(!el->HasDependency()) {
			ek->ComputeDestinationIndices();
			
			if(mineq != -1 || maxeq != -1)
			{
				FilterEquations(ek->fSourceIndex,ek->fDestinationIndex,mineq,maxeq);
			}
#ifdef LOG4CXX
			if(loggerel->isDebugEnabled())
			{
				std::stringstream sout;
				ek->fMat.Print("Element stiffness matrix",sout);
				ef->fMat.Print("Element right hand side", sout);
				LOGPZ_DEBUG(loggerel,sout.str())
			}
#endif
		} else {
			// the element has dependent nodes
			ek->ApplyConstraints();
			ef->ApplyConstraints();
			ek->ComputeDestinationIndices();
			if(mineq != -1 || maxeq != -1)
			{
				FilterEquations(ek->fSourceIndex,ek->fDestinationIndex,mineq,maxeq);
			}
		}
		
		
		// put the elementmatrices on the stack to be assembled (threadsafe)
		data->ComputedElementMatrix(iel,ek,ef);
		// compute the next element (this method is threadsafe)
		iel = data->NextElement();
	}
	return 0;
}

	// The function which will compute the assembly
void *TPZStructMatrix::ThreadData::ThreadAssembly(void *threaddata)
{
	ThreadData *data = (ThreadData *) threaddata;
	TPZCompMesh *cmesh = data->fMesh;
	int nel = cmesh->NElements();
	pthread_mutex_lock(&(data->fAccessElement));
	int nextel = data->fNextElement;
	int numprocessed = data->fProcessed.size();
	bool globalresult = true;
	while(nextel < nel || numprocessed)
	{
		std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > >::iterator itavail;
		std::set<int>::iterator itprocess;
		bool keeplooking = false;
		if(data->fSubmitted.size() && data->fProcessed.size())
		{
			itavail = data->fSubmitted.begin();
			itprocess = data->fProcessed.begin();
			if(itavail->first == *itprocess)
			{
				// make sure we come back to look for one more element
				keeplooking = true;
				// Get a hold of the data
				int iel = *itprocess;
				data->fProcessed.erase(itprocess);
				TPZAutoPointer<TPZElementMatrix> ek = itavail->second.first;
				TPZAutoPointer<TPZElementMatrix> ef = itavail->second.second;
				data->fSubmitted.erase(itavail);
#ifdef LOG4CXX
				std::stringstream sout;
				sout << "Assembling element " << iel;
				LOGPZ_DEBUG(logger,sout.str())
#endif
#ifdef CHECKCONSISTENCY
				//static TPZCheckConsistency stiffconsist("ElementStiff");
				stiffconsist.SetOverWrite(true);
				bool result;
				result = stiffconsist.CheckObject(ek->fMat);
				if(!result)
				{
					globalresult = false;
					std::stringstream sout;
					sout << "element " << iel << " computed differently";
					LOGPZ_ERROR(loggerCheck,sout.str())
				}
#endif
				// Release the mutex
				pthread_mutex_unlock(&data->fAccessElement);
				// Assemble the matrix
				if(!ek->HasDependency())
				{
					data->fGlobMatrix->AddKel(ek->fMat,ek->fSourceIndex,ek->fDestinationIndex);
					data->fGlobRhs->AddFel(ef->fMat,ek->fSourceIndex,ek->fDestinationIndex);				
				}
				else
				{
					data->fGlobMatrix->AddKel(ek->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
					data->fGlobRhs->AddFel(ef->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);				
				}
				// acquire the mutex
				pthread_mutex_lock(&data->fAccessElement);
			}
		}
		if(!keeplooking)
		{
			pthread_mutex_unlock(&data->fAccessElement);
			LOGPZ_DEBUG(logger,"Going to sleep within assembly")
			// wait for a signal
			sem_wait(data->fAssembly);
			LOGPZ_DEBUG(logger,"Waking up for assembly")
			pthread_mutex_lock(&data->fAccessElement);
		}
		nextel = data->fNextElement;
		numprocessed = data->fProcessed.size();

	}
//	std::cout << std::endl;
	{
		std::stringstream sout;
		sout << "nextel = " << nextel << " numprocessed = " << numprocessed << " submitted " << data->fSubmitted.size() << std::endl;
		sout << "The comparaison results are : consistency check " << globalresult;
		LOGPZ_DEBUG(loggerCheck,sout.str())
	}
	pthread_mutex_unlock(&data->fAccessElement);
	return 0;	
}		
	
int TPZStructMatrix::ThreadData::NextElement()
{
	pthread_mutex_lock(&fAccessElement);
	int iel;
	int nextel = fNextElement;
	TPZCompMesh *cmesh = fMesh;
	TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
	int nel = elementvec.NElements();
	for(iel=fNextElement; iel < nel; iel++)
	{
		TPZCompEl *el = elementvec[iel];
		if(!el) continue;
		if(fMaterialIds.size() == 0) break;
		TPZAutoPointer<TPZMaterial> mat = el->Material();
		TPZSubCompMesh *cmesh = dynamic_cast<TPZSubCompMesh *> (el);
		if(!mat) 
		{
			if(!cmesh)
			{
				continue;
			}
		}
		else 
		{
			int matid = mat->Id();
			if(fMaterialIds.find(matid) == fMaterialIds.end()) continue;
		}
		break;
	}
	fNextElement = iel+1;
	nextel = iel;
	if(iel<nel) fProcessed.insert(iel);
	pthread_mutex_unlock(&fAccessElement);
	return nextel;
}

	// put the computed element matrices in the map
void TPZStructMatrix::ThreadData::ComputedElementMatrix(int iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef)
{
	pthread_mutex_lock(&fAccessElement);
	std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > el(ek,ef);
	fSubmitted[iel] = el;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Putting the stiffness matrix of element " << iel << " on the stack, signaling the semaphore";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	sem_post(fAssembly);
	pthread_mutex_unlock(&fAccessElement);	
	
}

	/// Set the set of material ids which will be considered when assembling the system
	void TPZStructMatrix::SetMaterialIds(const std::set<int> &materialids)
	{
		fMaterialIds = materialids;
		{
			std::set<int>::const_iterator it;
			std::stringstream sout;
			sout << "setting input material ids ";
			for(it=materialids.begin(); it!= materialids.end(); it++)
			{
				sout << *it << " ";
			}
			LOGPZ_DEBUG(logger,sout.str())
		}
		if(!fMesh)
		{
			LOGPZ_WARN(logger,"SetMaterialIds called without mesh")
			return;
		}
		int iel;
		TPZAdmChunkVector<TPZCompEl*> &elvec = fMesh->ElementVec();
		int nel = elvec.NElements();
		for(iel=0; iel<nel; iel++)
		{
			TPZCompEl *cel = elvec[iel];
			if(!cel) continue;
			TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (cel);
			if(!subcmesh) continue;
			TPZAnalysis *anal = subcmesh->GetAnalysis();
			if(!anal)
			{
				LOGPZ_WARN(logger,"SetMaterialIds called for substructure without analysis object")
				subcmesh->SetAnalysis();
				anal=subcmesh->GetAnalysis();
			}
			TPZAutoPointer<TPZStructMatrix> str = anal->StructMatrix();
			if(!str)
			{
				LOGPZ_WARN(logger,"SetMaterialIds called for substructure without structural matrix")
				continue;
			}
			str->SetMaterialIds(materialids);
		}
	}
