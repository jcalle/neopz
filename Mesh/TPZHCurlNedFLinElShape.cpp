/**
* @file
* @brief Contains the implementation of the TPZHCurlNedFLinEl::Shape method.
*/

#include "TPZHCurlNedFLinEl.h"
#ifdef HCURL_HIERARCHICAL
#include "pzshapelinear.h"

using namespace pzshape;

void TPZHCurlNedFLinEl::CalcShape(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi,
                              TPZFMatrix<REAL> &curlPhiHat, TPZVec<int> &order, TPZVec<int> nShapeF){

    

    const int nCon = order.size();
    const int dim = TPZShapeLinear::Dimension;
    const int firstSide = TPZShapeLinear::NSides - TPZShapeLinear::NFaces - 1;

#ifdef PZDEBUG
    if (!(nCon == 1 && dim == 1 && firstSide == 0)) {
        DebugStop();
    }
#endif

    const int lastFuncPos = nShapeF[nCon - 1] - 1;

    phi.Resize(lastFuncPos + 1, dim);
    curlPhiHat.Resize(1, 1); // The curl wont be calculated in boundary
                             // elements for now.
    const int pOrder = order[0];
    int currentFuncPos = lastFuncPos;
    const REAL x = (1. + qsi[0])/2;
    switch (pOrder) {
    case 15:
        phi(currentFuncPos, 0) =
            (1163381400 *
             (1 - 210 * x + 10920 * pow(x, 2) -
              247520 * pow(x, 3) + 3063060 * pow(x, 4) -
              23279256 * pow(x, 5) + 116396280 * pow(x, 6) -
              399072960 * pow(x, 7) + 960269310 * pow(x, 8) -
              1636014380 * pow(x, 9) + 1963217256 * pow(x, 10) -
              1622493600 * pow(x, 11) + 878850700 * pow(x, 12) -
              280816200 * pow(x, 13) + 40116600 * pow(x, 14)));
        currentFuncPos--;
    case 14:
        phi(currentFuncPos, 0) =
            (280816200 *
             (-1 + 182 * x - 8190 * pow(x, 2) +
              160160 * pow(x, 3) - 1701700 * pow(x, 4) +
              11027016 * pow(x, 5) - 46558512 * pow(x, 6) +
              133024320 * pow(x, 7) - 261891630 * pow(x, 8) +
              355655300 * pow(x, 9) - 327202876 * pow(x, 10) +
              194699232 * pow(x, 11) - 67603900 * pow(x, 12) +
              10400600 * pow(x, 13)));
        currentFuncPos--;
    case 13:
        phi(currentFuncPos, 0) =

            (67603900 *
             (1 - 156 * x + 6006 * pow(x, 2) -
              100100 * pow(x, 3) + 900900 * pow(x, 4) -
              4900896 * pow(x, 5) + 17153136 * pow(x, 6) -
              39907296 * pow(x, 7) + 62355150 * pow(x, 8) -
              64664600 * pow(x, 9) + 42678636 * pow(x, 10) -
              16224936 * pow(x, 11) + 2704156 * pow(x, 12)));
        currentFuncPos--;
    case 12:
        phi(currentFuncPos, 0) =

            (16224936 * (-1 + 132 * x - 4290 * pow(x, 2) +
                         60060 * pow(x, 3) - 450450 * pow(x, 4) +
                         2018016 * pow(x, 5) - 5717712 * pow(x, 6) +
                         10501920 * pow(x, 7) - 12471030 * pow(x, 8) +
                         9237800 * pow(x, 9) - 3879876 * pow(x, 10) +
                         705432 * pow(x, 11)));
        currentFuncPos--;
    case 11:
        phi(currentFuncPos, 0) =

            (3879876 * (1 - 110 * x + 2970 * pow(x, 2) -
                        34320 * pow(x, 3) + 210210 * pow(x, 4) -
                        756756 * pow(x, 5) + 1681680 * pow(x, 6) -
                        2333760 * pow(x, 7) + 1969110 * pow(x, 8) -
                        923780 * pow(x, 9) + 184756 * pow(x, 10)));
        currentFuncPos--;
    case 10:
        phi(currentFuncPos, 0) =

            (923780 * (-1 + 90 * x - 1980 * pow(x, 2) +
                       18480 * pow(x, 3) - 90090 * pow(x, 4) +
                       252252 * pow(x, 5) - 420420 * pow(x, 6) +
                       411840 * pow(x, 7) - 218790 * pow(x, 8) +
                       48620 * pow(x, 9)));
        currentFuncPos--;
    case 9:
        phi(currentFuncPos, 0) =
            (218790 * (1 - 72 * x + 1260 * pow(x, 2) -
                       9240 * pow(x, 3) + 34650 * pow(x, 4) -
                       72072 * pow(x, 5) + 84084 * pow(x, 6) -
                       51480 * pow(x, 7) + 12870 * pow(x, 8)));
        currentFuncPos--;
    case 8:
        phi(currentFuncPos, 0) =
            (51480 *
             (-1 + 56 * x - 756 * pow(x, 2) + 4200 * pow(x, 3) -
              11550 * pow(x, 4) + 16632 * pow(x, 5) -
              12012 * pow(x, 6) + 3432 * pow(x, 7)));
        currentFuncPos--;
    case 7:
        phi(currentFuncPos, 0) =
            (12012 * (1 - 42 * x + 420 * pow(x, 2) -
                      1680 * pow(x, 3) + 3150 * pow(x, 4) -
                      2772 * pow(x, 5) + 924 * pow(x, 6)));
        currentFuncPos--;
    case 6:
        phi(currentFuncPos, 0) =
            (2772 *
             (-1 + 30 * x - 210 * pow(x, 2) + 560 * pow(x, 3) -
              630 * pow(x, 4) + 252 * pow(x, 5)));
        currentFuncPos--;
    case 5:
        phi(currentFuncPos, 0) =
            (630 * (1 - 20 * x + 90 * pow(x, 2) -
                    140 * pow(x, 3) + 70 * pow(x, 4)));
        currentFuncPos--;
    case 4:
        phi(currentFuncPos, 0) =
            (140 *
             (-1 + 12 * x - 30 * pow(x, 2) + 20 * pow(x, 3)));
        currentFuncPos--;
    case 3:
        phi(currentFuncPos, 0) = (30 * (1 - 6 * x + 6 * pow(x, 2)));
        currentFuncPos--;
    case 2:
        phi(currentFuncPos, 0) = (-6 + 12 * x);
        currentFuncPos--;
    case 1:
        phi(currentFuncPos, 0) = (1);
        break;
    default:
        DebugStop(); // polynomial order not implemented!
    }
}
#endif
