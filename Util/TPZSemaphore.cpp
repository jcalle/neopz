/** 
 * @file
 * @brief Contains the implementation of the methods to TPZSemaphore class.
 */

#include "TPZSemaphore.h"

#include <fstream>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.util.semaphore"));
#endif


TPZSemaphore::TPZSemaphore(int initCount)
{
	fCounter = initCount;
}

TPZSemaphore& TPZSemaphore::operator=(const TPZSemaphore &cp)
{
  fCounter = cp.fCounter;
  return *this;
}

TPZSemaphore::TPZSemaphore(const TPZSemaphore &cp){ *this = cp;}

void TPZSemaphore::Wait()
{
  std::unique_lock lck (fMutex);
	if (fCounter > 0)
	{
		fCounter--;
#ifdef LOG4CXX
		if(logger->isDebugEnabled())
		{
			std::stringstream sout;
#ifdef VC
			sout << "THREAD IN SEMAPHORE WAIT: " << pthread_self().x <<" " <<__LINE__ << std::endl;
#else
			sout << "THREAD IN SEMAPHORE WAIT: " << pthread_self() <<" " <<__LINE__ << std::endl;
#endif
			sout << "FCOUNTER VALUE : " << fCounter << std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return;
	}
	while (fCounter == 0) {
    fCond.wait(lck);
		if (fCounter > 0) {
			fCounter--;
#ifdef LOG4CXX
			if(logger->isDebugEnabled())
			{
				std::stringstream sout;
#ifdef VC
				sout << "THREAD IN SEMAPHORE AFTER WAIT: " << pthread_self().x <<" " << __LINE__ << std::endl;
#else
				sout << "THREAD IN SEMAPHORE AFTER WAIT: " << pthread_self() <<" " << __LINE__ << std::endl;
#endif
				sout << "FCOUNTER VALUE : " << fCounter << std::endl;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			return;
		}
	}
	
}

void TPZSemaphore::Post()
{
	std::unique_lock lck(fMutex);
        
	fCounter++;
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
#ifdef VC
		sout << "THREAD IN SEMAPHORE POST: " << pthread_self().x <<" " << __LINE__ << std::endl;
#else
		sout << "THREAD IN SEMAPHORE POST: " << pthread_self() <<" " << __LINE__ << std::endl;
#endif
		sout << "FCOUNTER VALUE : " << fCounter << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  lck.unlock();
  fCond.notify_one();
	return;
}
