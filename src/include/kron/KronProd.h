#ifndef KRONPROD_H
#define KRONPROD_H

#include <vector>
#include <armadillo>
#include "include/spin/SpinInteractionDefine.h"

using namespace std;
using namespace arma;

/// \defgroup KronProd Kron
///@{

////////////////////////////////////////////////////////////////////////////////
//{{{ KronProd
class KronProd
{
public:
    KronProd();
    KronProd(DIM_LIST dim_list);
    ~KronProd();

    void   fill(INDICES idx, MULTIPLIER coeff, TERM mat);
    cx_mat full();
    DIM_LIST getDimList(){return _dim_list;};

    friend KronProd Flat(const KronProd& kp);
    friend KronProd Sharp(const KronProd& kp);
    friend KronProd CircleC(const KronProd& kp);

    friend ostream&  operator << (ostream& outs, const KronProd& kp);
protected:
private:
    DIM_LIST   _dim_list;
    INDICES    _spin_index;
    MULTIPLIER _coeff;
    TERM       _mat;
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SumKronProd
class SumKronProd
{
public:
    SumKronProd();
    SumKronProd(const vector<KronProd>& kp_lst);
    ~SumKronProd();

    cx_mat full();
    vector<KronProd> getKronProdList(){return _kron_prod_list;};
    DIM_LIST getDimList(){return _dim_list;};
    
    friend SumKronProd Flat(const SumKronProd& skp);
    friend SumKronProd Sharp(const SumKronProd& skp);
    friend SumKronProd CircleC(const SumKronProd& skp);
    friend SumKronProd& operator + (SumKronProd& sum, const SumKronProd skp);
    friend ostream&  operator << (ostream& outs, SumKronProd& skp);
protected:
    DIM_LIST _dim_list;
    vector<KronProd> _kron_prod_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

SumKronProd Flat(const SumKronProd& skp);
SumKronProd Sharp(const SumKronProd& skp);
SumKronProd CircleC(const SumKronProd& skp);

typedef SumKronProd(FUNC)(const SumKronProd&);

extern FUNC* FlatOperation;
extern FUNC* SharpOperation;
extern FUNC* CircleCOperation;
/// @}
#endif

