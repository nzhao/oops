#include <armadillo>
#include "include/kron/KronProd.h"
#include "include/easylogging++.h"


////////////////////////////////////////////////////////////////////////////////
//{{{ MatrixOperationFunctions
//using namespace std::placeholders;

cx_mat Flat(const cx_mat& m)
{
    mat id; id.eye(size(m));
    return  kron(id, m);
}

cx_mat Sharp(const cx_mat& m)
{
    mat id; id.eye(size(m));
    return  kron(m.st(), id);
}

cx_mat CircleC(const cx_mat& m)
{
    mat id; id.eye(size(m));
    return   kron(id,m) - kron(conj(m),id);
}
MatExpanFunc* FLAT = &Flat;
MatExpanFunc* SHARP = &Sharp;
MatExpanFunc* CIRCLEC = &CircleC;
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ KronProd

KronProd::KronProd()
{ //LOG(INFO) << "Default constructor: KronProd";
}

KronProd::~KronProd()
{ //LOG(INFO) << "Default destructor: KronProd";
}

KronProd::KronProd(DIM_LIST dim_list)
{
    _dim_list = dim_list;
    _kron_num = dim_list.size();

    _dim  = 1;
    for(int i=0; i<_dim_list.size(); ++i)
        _dim *= _dim_list[i];
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

cx_vec KronProd::vecterize()
{
    vector<cx_mat> all_mat;
    for(int i=0; i<_dim_list.size(); ++i)
        all_mat.push_back( eye<cx_mat>(_dim_list[i], _dim_list[i]) );

    for(int i=0; i<_spin_index.size(); ++i)
        all_mat[ _spin_index[i] ] = _mat[i];

    cx_vec res=vectorise(all_mat[0]);
    for(int i=1; i<all_mat.size(); ++i)
        res=kron(res, vectorise(all_mat[i]) );

    return _coeff*res;
}

KronProd Expand(const KronProd& kp, MatExpanFunc* expan_func)
{
    DIM_LIST new_dim;
    //for (auto d : kp._dim_list)
    //    new_dim.push_back( d*d );
    for (int i=0; i<kp._dim_list.size(); ++i)
        new_dim.push_back( kp._dim_list[i]*kp._dim_list[i] );
    KronProd res(new_dim);

    TERM new_mat;
    //for(cx_mat A : kp._mat)
    //{
    //    mat id; id.eye(size(A));
    //    new_mat.push_back( expan_func(A) );
    //}
    for(int i=0; i<kp._mat.size(); ++i)
    {
        cx_mat A( kp._mat[i] );
        mat id; id.eye(size(A));
        new_mat.push_back( expan_func(A) );
    }
    res.fill(kp._spin_index, kp._coeff, new_mat);
    return res;
}

ostream&  operator << (ostream& outs, const KronProd& kp)
{
    outs << "=======================================================================================================================================================================================" << endl;
    outs << "DIM_LIST= ";
    //for(auto dim:kp._dim_list)
    //    outs << dim << ", ";
    for(int i=0; i<kp._dim_list.size(); ++i)
        outs << kp._dim_list[i] << ", ";
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
{ //LOG(INFO) << "Default constructor: SumKronProd.";
}

SumKronProd::SumKronProd(const vector<KronProd>& kp_lst)
{ //LOG(INFO) << "Constructor: SumKronProd with KronProd list.";
    _kron_prod_list = kp_lst;
    if( !_kron_prod_list.empty() )
    {
        _dim_list = _kron_prod_list[0].getDimList();
        _dim = _kron_prod_list[0].getDim();
        _kron_num = _kron_prod_list[0].getKronNum();
    }
}
SumKronProd::~SumKronProd()
{ //LOG(INFO) << "Default destructor: SumKronProd.";
}
cx_mat SumKronProd::full()
{
    cx_mat res = _kron_prod_list[0].full();
    for(int i=1; i<_kron_prod_list.size(); ++i)
        res = res + _kron_prod_list[i].full();
    return res;
}

cx_vec SumKronProd::vecterize()
{
    cx_vec res = _kron_prod_list[0].vecterize();
    for(int i=1; i<_kron_prod_list.size(); ++i)
        res = res + _kron_prod_list[i].vecterize();
    return res;
}

SumKronProd Expand(const SumKronProd& skp, MatExpanFunc* expan_func)
{
    vector<KronProd> new_kp_list;
    new_kp_list.reserve( skp._kron_prod_list.size() );
    //for( auto kp : skp._kron_prod_list)
    //    new_kp_list.push_back( Expand(kp, expan_func) );
    for(int i=0; i<skp._kron_prod_list.size(); ++i)
        new_kp_list.push_back( Expand(skp._kron_prod_list[i], expan_func) );

    SumKronProd res(new_kp_list);
    return res;
}

vector<MULTIPLIER> SumKronProd::getCoeffList() const
{
    vector<MULTIPLIER> res;
    for(int i=0; i<_kron_prod_list.size(); ++i)
        res.push_back( _kron_prod_list[i].getCoeff() );
    return res;
}

vector<INDICES> SumKronProd::getIndicesList() const
{
    vector<INDICES> res;
    for(int i=0; i<_kron_prod_list.size(); ++i)
        res.push_back( _kron_prod_list[i].getIndices() );
    return res;
}

vector<TERM> SumKronProd::getTermList() const
{
    vector<TERM> res;
    for(int i=0; i<_kron_prod_list.size(); ++i)
        res.push_back( _kron_prod_list[i].getTermMat() );
    return res;
}
vector<int> SumKronProd::getMatNumList() const
{
    vector<int> res;
    for(int i=0; i<_kron_prod_list.size(); ++i)
        res.push_back( _kron_prod_list[i].getMatNum() );
    return res;
}

SumKronProd& SumKronProd::scale(double factor)
{
    for(vector<KronProd>::iterator it=_kron_prod_list.begin(); it != _kron_prod_list.end(); ++it)
        it->scale(factor);
    return *this;
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
    //for(auto kp:skp._kron_prod_list)
    //{
    //    outs << kp << endl;
    //    outs << "KronProd.full() = " <<  endl;
    //    outs << kp.full() << endl;
    //    outs << endl;
    //}
    for(int i=0; i<skp._kron_prod_list.size(); ++i)
    {
        outs << skp._kron_prod_list[i] << endl;
        outs << "KronProd.full() = " <<  endl;
        outs << skp._kron_prod_list[i].full() << endl;
        outs << endl;
    }
    outs << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    outs << "Sum Matrix = " << endl;
    outs << skp.full() << endl;
    return outs;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
