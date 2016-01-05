#include <iostream>
#include <armadillo>
#include "include/spin/SpinClusterAlgorithm.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////////////////////////
//{{{ cClusterIndex
cClusterIndex::cClusterIndex()
{// cout << "cClusterIndex default constructor" << endl;
}

cClusterIndex::cClusterIndex(const uvec& idx)
{
    _index = idx;
    sort(_index.begin(), _index.end()); }

cClusterIndex::~cClusterIndex()
{ //cout << "cClusterIndex default destructor" << endl;
}

mat cClusterIndex::get_array(size_t nspin)
{
    mat idx_array=zeros(1, nspin);
    int nnz = _index.size();

    for(int i=0; i<nnz; ++i)
        idx_array[_index[i]]=1;
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
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinGrouping
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
    int nClst=_cluster_index_list[order].size();
    mat res=zeros(nClst, _nspin);

    int i=0;
    for(cClusterIndex vIdx:_cluster_index_list[order])
    {
        res.row(i) = vIdx.get_array(_nspin);
        i++;
    }
    return conv_to<sp_mat>::from(res);
}

void cSpinGrouping::subgraph2index(const sp_mat& subgraph)
{
    for(int i=0; i<subgraph.n_rows; ++i)
    {
        cout << "\r" <<  i << "/" << subgraph.n_rows << "subgraphs are inserted. \t";
        mat r(subgraph.row(i));  uvec nz_r = find(r);  int order = nz_r.size()-1;
        cClusterIndex cIdx( nz_r );
        _cluster_index_list[ order ].insert(cIdx);
        //_cluster_index_list[ order ].push_back(cIdx);
    }
    cout <<endl;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinDepthFirstPathTracing
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
    mat new_subgraph=zeros( sum(nonzeros(neighbor)), _nspin);

    int nGen=0;
    for(int i=0; i<subgraph.n_rows; ++i)
    {
        cout << "\r"<< i+1 <<  "/" << subgraph.n_rows
             << " clusters are generated of spin order= " << subgraph_order+1 << "\t";

        mat parent_row(subgraph.row(i));
        mat r(neighbor.row(i));    uvec candidate = find(r);

        int row_idx=0;
        for(int j=0; j<candidate.size(); ++j)
        {
            if(parent_row(candidate(j))==0)
            {
                new_subgraph.row(nGen)=parent_row;
                new_subgraph(nGen, candidate(j))=1;
                nGen++;
            }
        }
        cout.flush();
    }
    cout << endl;

    sp_mat res_mat=conv_to<sp_mat>::from( new_subgraph.rows(0, nGen-1) );

    return res_mat;
    //return remove_repeat(res_mat, subgraph_order);
}

sp_mat cDepthFirstPathTracing::remove_repeat(sp_mat subgraph, int subgraph_order)
{
    cout << "removing repeated clusters... " << endl;
    int nGen=subgraph.n_rows;
    sp_mat res_mat2 = subgraph*subgraph.t();
    for(int i=0; i<nGen; ++i)
        res_mat2(i, i)=0;
    cout << "mat prod finished." << endl;

    set<int> to_remove;
    for(sp_mat::const_iterator it = res_mat2.begin(); it != res_mat2.end(); ++it)
    {
        if( (*it) == subgraph_order+1 )
        {
            to_remove.insert( it.row() > it.col() ? it.row() : it.col() );
        }
    }

    vec idx=zeros(subgraph.n_rows);
    int count =0;
    for(int i=0; i<nGen; ++i)
    {
        if( to_remove.find(i) == to_remove.end() )
        {
            idx(count)=i;
            count++;
        }
    }
    cout << "finished." << endl;

    uvec x=conv_to<uvec>::from(idx.rows(0, count-1));
    mat subgraph_full=conv_to<mat>::from( subgraph );
    sp_mat res_mat_sel=conv_to<sp_mat>::from( subgraph_full.rows(x) );
    return res_mat_sel;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
