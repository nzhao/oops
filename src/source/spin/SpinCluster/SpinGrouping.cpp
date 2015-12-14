#include <iostream>
#include <armadillo>
#include "include/spin/SpinGrouping.h"

using namespace arma;


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
// cClusterIndex
cClusterIndex::cClusterIndex()
{
    cout << "cClusterIndex default constructor" << endl;
}

cClusterIndex::cClusterIndex(const uvec& idx)
{
    _index = idx;
    sort(_index.begin(), _index.end());
    cout << "class properties" << &_index << "\t local variable" << &idx << endl;
}

cClusterIndex::~cClusterIndex()
{
    cout << "cClusterIndex default destructor" << endl;
}

sp_mat cClusterIndex::get_array(size_t nspin)
{
    int nnz = _index.size();
    vec values=ones<vec>(nnz);
    umat locations=zeros<umat>(2,nnz);

    for(int i=0; i<nnz; ++i)
        locations(1, i)=_index[i];

    sp_mat idx_array(locations, values, 1, nspin);
    return idx_array;
}

bool operator == (const cClusterIndex& idx1, const cClusterIndex& idx2)
{
    if(idx1._index.size() == idx2._index.size())
        for(int i=0; i<idx1._index.size(); ++i)
        { if(idx1._index[i] != idx2._index[i]) return 0; }
    else
        return 0;
    return 1;
}


bool operator < (const cClusterIndex& idx1, const cClusterIndex& idx2)
{
    int sz1=idx1._index.size(); int sz2=idx2._index.size();
    if( sz1 == sz2 )
        for(int i=0; i<sz1; ++i)
        { if(idx1._index[i] != idx2._index[i]) return idx1._index[i] < idx2._index[i]; }
    else
        return (sz1 < sz2);
    return 0;// idx1 equals to idx2
}

ostream&  operator << (ostream& outs, const cClusterIndex& idx)
{
    for(auto i: idx._index)
        outs << i << ", ";
    return outs;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinGrouping
cSpinGrouping::cSpinGrouping()
{
    cout << "cSpinGrouping default constructor is called." << endl;
    _cluster_index_list = {};
}
cSpinGrouping::cSpinGrouping(umat connection_matrix)
{
    _connection_matrix=connection_matrix;
    cout << "cSpinGrouping constructor is called." << endl;
}

cSpinGrouping::~cSpinGrouping()
{
    cout << "cSpinGrouping destructor is called." << endl;
}

arma::sp_mat cSpinGrouping::get_cluster_mat(int order)
{
    int nspin=_connection_matrix.size();
    sp_mat res={};
    for(cClusterIndex vIdx:_cluster_index_list[order])
        res=join_cols(res, vIdx.get_array(nspin));
    return res;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinDepthFirstPathTracing
cDepthFirstPathTracing::cDepthFirstPathTracing()
{

    cout << "need a connection matrix." << endl;
}

cDepthFirstPathTracing::cDepthFirstPathTracing(umat  connection_matrix)
{
    cout << "cSpinDepthFirstPathTracing constructor is called." << endl;
    _connection_matrix=connection_matrix;
}

cDepthFirstPathTracing::~cDepthFirstPathTracing()
{
    cout << "cSpinDepthFirstPathTracing destructor is called." << endl;
}

CLST cDepthFirstPathTracing::generate()
{
    mat subgraph = eye<mat>(size(_connection_matrix));
    sp_mat sp_subgraph(subgraph);
    cout << sp_subgraph << endl;
    cout << "list length= " <<  _cluster_index_list.size() << endl;

    sp_subgraph(0,5)=1;

    size_t iii;
    sp_mat X=sp_subgraph.row(0);
    sp_mat::const_iterator start = X.begin(); sp_mat::const_iterator end = X.end();
    cout << "nonzero idx:" << endl;
    for( sp_mat::const_iterator it = start; it != end; ++it)
    {       //cout << sp_subgraph.row(0) << endl;
        iii=it.col();
        cout << iii << endl;
    }
//    uvec qq = find(sp_subgraph.row(0));
//    cout << "qq" << qq << endl;
    return _cluster_index_list;
}
