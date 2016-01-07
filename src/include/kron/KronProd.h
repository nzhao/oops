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

    KronProd Flat();
    KronProd Sharp();
    KronProd FlatSharp();
    KronProd CircleC();

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
    friend SumKronProd& operator + (SumKronProd& sum, const SumKronProd skp);
    friend ostream&  operator << (ostream& outs, SumKronProd& skp);
protected:
    vector<KronProd> _kron_prod_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
#endif

