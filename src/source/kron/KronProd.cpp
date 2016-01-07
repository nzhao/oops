#include <armadillo>
#include "include/kron/KronProd.h"
#include "include/easylogging++.h"

using namespace std::placeholders;

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

KronProd Flat(const KronProd& kp)
{
    DIM_LIST new_dim;
    transform(kp._dim_list.cbegin(), kp._dim_list.cend(),
              new_dim.begin(), bind(multiplies<int>(), _1, _1));
    KronProd res(new_dim);

    TERM new_mat;
    for(cx_mat A : kp._mat)
    {
        mat id; id.eye(size(A));
        new_mat.push_back( kron(id, A) );
    }

    res.fill(kp._spin_index, kp._coeff, new_mat);
    return res;
}

KronProd Sharp(const KronProd& kp)
{
    DIM_LIST new_dim;
    transform(kp._dim_list.cbegin(), kp._dim_list.cend(),
              new_dim.begin(), bind(multiplies<int>(), _1, _1));
    KronProd res(new_dim);

    TERM new_mat;
    for(cx_mat B : kp._mat)
    {
        mat id; id.eye(size(B));
        new_mat.push_back( kron(B.st(), id) );
    }

    res.fill(kp._spin_index, kp._coeff, new_mat);
    return res;
}

KronProd CircleC(const KronProd& kp)
{
    DIM_LIST new_dim;
    transform(kp._dim_list.cbegin(), kp._dim_list.cend(),
              new_dim.begin(), bind(multiplies<int>(), _1, _1));
    KronProd res(new_dim);

    TERM new_mat;
    for(cx_mat H : kp._mat)
    {
        mat id; id.eye(size(H));
        new_mat.push_back( kron(id,H) - kron(conj(H),id) );
    }

    res.fill(kp._spin_index, kp._coeff, new_mat);
    return res;
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

SumKronProd Flat(const SumKronProd& skp)
{
    vector<KronProd> new_kp_list;
    transform(skp._kron_prod_list.cbegin(), skp._kron_prod_list.cend(),
              new_kp_list.begin(),
              [](KronProd kp) { return Flat(kp);});

    SumKronProd res(new_kp_list);
    return res;
}

SumKronProd Sharp(const SumKronProd& skp)
{
    vector<KronProd> new_kp_list;
    transform(skp._kron_prod_list.cbegin(), skp._kron_prod_list.cend(),
              new_kp_list.begin(),
              [](KronProd kp) { return Sharp(kp);});

    SumKronProd res(new_kp_list);
    return res;
}

SumKronProd CircleC(const SumKronProd& skp)
{
    vector<KronProd> new_kp_list;
    transform(skp._kron_prod_list.cbegin(), skp._kron_prod_list.cend(),
              new_kp_list.begin(),
              [](KronProd kp) { return CircleC(kp);});

    SumKronProd res(new_kp_list);
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
