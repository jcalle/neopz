/* 
 * File:   TPZConstrainedNormalRandom.h
 * Author: quinelato
 *
 * Created on 5 de Dezembro de 2017, 15:28
 */

#ifndef TPZCONSTRAINEDNORMALRANDOM_H
#define TPZCONSTRAINEDNORMALRANDOM_H

#include "TPZConstrainedRandom.h"
#include "TPZNormalRandom.h"

template <typename TVar>
class TPZConstrainedNormalRandom : public TPZConstrainedRandom<TVar>, public TPZNormalRandom<TVar> {
public:
    TPZConstrainedNormalRandom(TVar begin, TVar end, TVar mean, TVar stdev): TPZConstrainedRandom(begin, end), TPZNormalRandom(mean, stdev) {
    }
    TPZConstrainedNormalRandom(const TPZConstrainedNormalRandom<TVar>& orig): TPZConstrainedRandom(orig), TPZNormalRandom(orig) {
    }
    virtual TPZRandom<TVar> *clone(){
        return new TPZConstrainedNormalRandom<TVar>(*this);
    }
    virtual TVar next();
    TVar pdf(TVar x);
    virtual ~TPZConstrainedNormalRandom(){
    }

};

template <typename TVar>
TVar TPZConstrainedNormalRandom<TVar>::next() {
    TVar value;
    do {
        value = TPZNormalRandom<TVar>::next();
    } while (value <= begin || value >= end);
    return value;
}

template <typename TVar>
TVar TPZConstrainedNormalRandom<TVar>::pdf(TVar x) {
    if (x<begin || x > end) return 0;
    TVar normal_pdf = TPZNormalRandom<TVar>::pdf(x);
    TVar area = TPZNormalRandom<TVar>::cdf(end)-TPZNormalRandom<TVar>::cdf(begin);
    return normal_pdf / area;
}

#endif /* TPZCONSTRAINEDNORMALRANDOM_H */

