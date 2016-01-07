#include <armadillo>
#include "include/kron/KronProd.h"

#include "include/easylogging++.h"
////////////////////////////////////////////////////////////////////////////////
//{{{ KronProd

KronProd::KronProd()
{ LOG(INFO) << "Default constructor: KronProd";}

KronProd::~KronProd()
{ LOG(INFO) << "Default destructor: KronProd";}

KronProd::KronProd(DIM_LIST dim_list)
{
    _dim_list = dim_list;
}

void KronProd::fill(INDICES idx, MULTIPLIER coeff, TERM mat)
{
    _spin_index = idx;
    _coeff      = coeff;
    _mat        = mat;
}

cx_mat KronProd::full()
{
    vector<cx_mat> all_mat;
    for(int i=0; i<_dim_list.size(); ++i)
        all_mat.push_back( eye<cx_mat>(_dim_list[i], _dim_list[i]) );

    for(int i=0; i<_spin_index.size(); ++i)
        all_mat[ _spin_index[i] ] = _mat[i];

    cx_mat res=all_mat[0];
    for(int i=1; i<all_mat.size(); ++i)
        res=kron(res, all_mat[i]);

    return _coeff*res;
}

KronProd& KronProd::Flat()
{
}

ostream&  operator << (ostream& outs, const KronProd& kp)
{
    outs << "=======================================================================================================================================================================================" << endl;
    outs << "DIM_LIST= ";
    for(auto dim:kp._dim_list)
        outs << dim << ", ";
    outs << endl;

    outs << "COEFF= " << kp._coeff << endl;

    for(int i=0; i<kp._spin_index.size(); ++i)
    {
        outs << "SPIN[" <<kp._spin_index[i] << "] = " << endl;
        outs << kp._mat[i] << endl;
    }

    return outs;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SumKronProd
SumKronProd::SumKronProd()
{
    cout << "Default SumKronProd constructor" << endl;
}

SumKronProd::SumKronProd(const vector<KronProd>& kp_lst)
{
    cout << "SumKronProd constructor with kp_lst" << endl;
    _kron_prod_list = kp_lst;
}
SumKronProd::~SumKronProd()
{
    cout << "Destructor of SumKronProd" << endl;
}
cx_mat SumKronProd::full()
{
    cx_mat res = _kron_prod_list[0].full();
    for(int i=1; i<_kron_prod_list.size(); ++i)
        res = res + _kron_prod_list[i].full();
    return res;
}
SumKronProd& operator + (SumKronProd& sum, const SumKronProd skp)
{
    vector<KronProd> A = sum._kron_prod_list;
    vector<KronProd> B = skp._kron_prod_list;
    vector<KronProd> AB;
    
    AB.reserve( A.size() + B.size() ); 
    AB.insert( AB.end(), A.begin(), A.end() );
    AB.insert( AB.end(), B.begin(), B.end() );
    sum._kron_prod_list = AB;

    return sum;
}
ostream&  operator << (ostream& outs, SumKronProd& skp)
{
    for(auto kp:skp._kron_prod_list)
    {
        outs << kp << endl;
        outs << "KronProd.full() = " <<  endl;
        outs << kp.full() << endl;
        outs << endl;
    }
    outs << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    outs << "Sum Matrix = " << endl;
    outs << skp.full() << endl;
    return outs;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
