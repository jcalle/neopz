#ifndef TPBrSimulationData_H
#define TPBrSimulationData_H

#include "TPBrStrainStressDataBase.h"
#include "TPZSandlerDimaggio.h"

class TPBrSimulationData : public TPBrStrainStressDataBase
{

public:
    TPBrSimulationData();
    TPBrSimulationData(int startidx, int endidx, int medidx);
    
    virtual ~TPBrSimulationData();
    
    TPBrSimulationData(const TPBrSimulationData &copy) : TPBrStrainStressDataBase(copy), fstart_idx(copy.fstart_idx),
        fend_idx(copy.fend_idx), fMedicao_idx(copy.fMedicao_idx), fSandler(copy.fSandler)
    {
        
    }
    
    TPBrSimulationData &operator=(const TPBrSimulationData &copy)
    {
        TPBrStrainStressDataBase::operator=(copy);
        fstart_idx = copy.fstart_idx;
        fend_idx = copy.fend_idx;
        fMedicao_idx = copy.fMedicao_idx;
        fSandler = copy.fSandler;
        return *this;
    }

    inline void Set_start_idx(int startidx) {
        fstart_idx = startidx;
    }
    inline void Set_end_idx(int endidx) {
        fend_idx = endidx;
    }
    inline int Get_start_idx(){
        return fstart_idx;
    }
    inline int Get_end_idx(){
        return fend_idx;
    }
    inline void Set_medicao_idx(int medidx){
        fMedicao_idx = medidx;
    }
    inline int Get_medicao_idx(){
        return fMedicao_idx;
    }
    inline void SetUpSandlerDimaggio(REAL poisson, REAL E, REAL A, REAL B, REAL C, REAL R, REAL D, REAL W)
    {
        fSandler.SetUp(poisson, E, A, B, C, R, D, W);
    }
    
    inline void SetUpSandlerDimaggio(const TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &copy)
    {
        fSandler = copy;
    }


protected:
    /// primeiro indice da medicao
    int fstart_idx;
    /// ultimo indice da medicao
    int fend_idx;
    /// Indice da medicao que deu origem a esta simulacao
    int fMedicao_idx;
    /// modelo plastico que gerou esta simulacao
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> fSandler;

};

#endif // TPBrSimulationData_H