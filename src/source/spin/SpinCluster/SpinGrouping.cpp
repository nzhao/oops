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

cClusterIndex::cClusterIndex(const vector<int>& idx)
{
    _index = idx;
    sort(_index.begin(), _index.end());
    cout << "class properties" << &_index << "\t local variable" << &idx << endl;
}

cClusterIndex::~cClusterIndex()
{
    cout << "cClusterIndex default destructor" << endl;
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
    return 1;
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

CLST_IDX_LIST cDepthFirstPathTracing::generate()
{
    mat subgraph = eye<mat>(size(_connection_matrix));
    cout << subgraph << endl;
    cout << "list length= " <<  _cluster_index_list.size() << endl;

    subgraph(0,5)=1;
    cout << subgraph << endl;
    uvec qq = find(subgraph.row(0));
    cout << "qq" << qq << endl;
    return _cluster_index_list;
}
