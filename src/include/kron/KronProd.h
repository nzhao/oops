#ifndef KRONPROD_H
#define KRONPROD_H

#include <vector>
#include <armadillo>
#include "include/kron/KronProd.h"
#include "include/spin/SpinInteractionDefine.h"

using namespace std;
using namespace arma;

/// \defgroup KronProd Kron
///@{

////////////////////////////////////////////////////////////////////////////////
//{{{ Matrix Operation Functions
cx_mat Flat(const cx_mat& m);
cx_mat Sharp(const cx_mat& m);
cx_mat CircleC(const cx_mat& m);
typedef cx_mat(MatExpanFunc) (const cx_mat&);

extern MatExpanFunc* FLAT;
extern MatExpanFunc* SHARP;
extern MatExpanFunc* CIRCLEC;
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ KronProd
class KronProd
{
public:
    KronProd();
    KronProd(DIM_LIST dim_list);
    ~KronProd();

    cx_mat      full();
    cx_vec      vecterize();
    void        fill(INDICES idx, MULTIPLIER coeff, TERM mat);
    KronProd&   scale(double factor) { _coeff *= factor; return *this;};
    DIM_LIST    getDimList(){return _dim_list;};
    int         getDim() const {return _dim;};
    MULTIPLIER  getCoeff() const {return _coeff;};
    INDICES     getIndices() const {return _spin_index;};
    TERM        getTermMat() const {return _mat;};
    int         getKronNum() const {return _kron_num;};
    int         getMatNum() const {return _mat.size();};

    friend KronProd Expand(const KronProd& kp, MatExpanFunc * exppan_func);
    friend ostream&  operator << (ostream& outs, const KronProd& kp);
protected:
private:
    int        _kron_num;
    DIM_LIST   _dim_list;
    INDICES    _spin_index;
    MULTIPLIER _coeff;
    TERM       _mat;
    int        _dim;
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
    cx_vec vecterize();
    vector<KronProd> getKronProdList(){return _kron_prod_list;};
    DIM_LIST getDimList(){return _dim_list;};
    int      getDim() const {return _dim;};
    int      getKronProdSize() const {return _kron_prod_list.size();};

    SumKronProd&  scale(double factor);
    void append(KronProd kp) {_kron_prod_list.push_back(kp);};
    int         getKronNum() const {return _kron_num;};

    vector<MULTIPLIER> getCoeffList() const;
    vector<INDICES> getIndicesList() const;
    vector<TERM>    getTermList() const;
    vector<int>     getMatNumList() const;
    
    friend SumKronProd Expand(const SumKronProd& skp, MatExpanFunc * exppan_func);
    friend SumKronProd& operator + (SumKronProd& sum, const SumKronProd skp);
    friend ostream&  operator << (ostream& outs, SumKronProd& skp);
protected:
    int      _dim;
    int      _kron_num;
    DIM_LIST _dim_list;
    vector<KronProd> _kron_prod_list;
    //vector<MULTIPLIER> _coeff_list;
    //vector<INDICES>    _indices_list;
    //vector<TERM>       _term_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
#endif

