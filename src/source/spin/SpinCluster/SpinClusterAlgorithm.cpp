#include <iostream>
#include <armadillo>
#include "include/spin/SpinClusterAlgorithm.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
// cClusterIndex
cClusterIndex::cClusterIndex()
{// cout << "cClusterIndex default constructor" << endl;
}

cClusterIndex::cClusterIndex(const uvec& idx)
{
    _index = idx;
    sort(_index.begin(), _index.end());
}

cClusterIndex::~cClusterIndex()
{ //cout << "cClusterIndex default destructor" << endl;
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
    for(auto it=idx._index.begin(); it !=idx._index.end(); ++it)
    {
        outs << *it;
        if(next(it) != idx._index.end())
            outs << ", ";
    }
    return outs;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinGrouping
cSpinGrouping::cSpinGrouping()
{
    _cluster_index_list = CLST_IDX_LIST(MAX_CLUSTER_ORDER);
}
cSpinGrouping::cSpinGrouping(const sp_mat& connection_matrix)
{// cout << "cSpinGrouping constructor is called." << endl;
    _connection_matrix=connection_matrix;
}

cSpinGrouping::~cSpinGrouping()
{// cout << "cSpinGrouping destructor is called." << endl;
}

sp_mat cSpinGrouping::index2subgraph(int order)
{
    sp_mat res={};
    for(cClusterIndex vIdx:_cluster_index_list[order])
        res=join_cols(res, vIdx.get_array(_nspin));
    return res;
}

void cSpinGrouping::subgraph2index(const sp_mat& subgraph)
{
    for(int i=0; i<subgraph.n_rows; ++i)
    {
        mat r(subgraph.row(i));  uvec nz_r = find(r);  int order = nz_r.size()-1;
        cClusterIndex cIdx( nz_r );
        _cluster_index_list[ order ].insert(cIdx);
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinDepthFirstPathTracing
cDepthFirstPathTracing::cDepthFirstPathTracing()
{
    cout << "need a connection matrix." << endl;
}

cDepthFirstPathTracing::cDepthFirstPathTracing(const sp_mat&  connection_matrix, size_t maxOrder)
{
    _max_order = maxOrder;
    _nspin     = connection_matrix.n_cols;
    _connection_matrix=connection_matrix;
    subgraph2index( speye(_nspin, _nspin) );
}

cDepthFirstPathTracing::~cDepthFirstPathTracing()
{// cout << "cSpinDepthFirstPathTracing destructor is called." << endl;
}

void cDepthFirstPathTracing::generate()
{
    sp_mat subgraph = index2subgraph(0);

    for( int i = 1; i < _max_order; ++i)
    {
        sp_mat neighbor = subgraph*_connection_matrix;
        sp_mat new_subgraph = subgraph_growth(subgraph, neighbor, i);
        subgraph2index(new_subgraph);
        subgraph= index2subgraph(i);
    }
}

sp_mat cDepthFirstPathTracing::subgraph_growth(const sp_mat& subgraph, const sp_mat& neighbor, int subgraph_order)
{
    sp_mat new_subgraph={};
    int n_child = 0;

    for(int i=0; i<subgraph.n_rows; ++i)
    {
        sp_mat parent_row(subgraph.row(i));
        mat r(neighbor.row(i));    uvec candidate = find(r);

        //cout << "p: " << conv_to<mat>::from(parent_row) << endl;
        for(int j=0; j<candidate.size(); ++j)
        {
            sp_mat temp_row=parent_row;
            //cout << candidate(j) << ": ";
            if( parent_row(candidate(j))==0 )
            {
                temp_row( candidate(j) ) = 1;
                new_subgraph=join_cols(new_subgraph, temp_row);
                n_child++;
                //cout << conv_to<mat>::from(temp_row);
            }
            //cout << endl;
        }
    }

    //cout << "new_subgraph:" << endl;
    //cout << conv_to<mat>::from(new_subgraph) << endl;
    return new_subgraph;
}
