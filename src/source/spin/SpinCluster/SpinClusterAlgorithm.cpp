#include <iostream>
#include <armadillo>
#include <iomanip> 
#include "include/spin/SpinClusterAlgorithm.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////////////////////////
//{{{ cClusterIndex
cClusterIndex::cClusterIndex()
{ //LOG(INFO) << "Default constructor of cClusterIndex.";
}

cClusterIndex::cClusterIndex(const uvec& idx)
{
    _index = idx;
    sort(_index.begin(), _index.end());
    _spin_num = idx.n_elem;
}

cClusterIndex::~cClusterIndex()
{ //LOG(INFO) << "Default destructor of cClusterIndex.";
}

mat cClusterIndex::get_array(size_t nspin)
{
    mat idx_array=zeros(1, nspin);
    size_t nnz = _index.size();

    for(int i=0; i<nnz; ++i)
        idx_array[_index[i]]=1;
    return idx_array;
}
set< CluserPostion > cClusterIndex::getSubClstPos() const
{
    set< CluserPostion > res;
    size_t order = getOrder();
    if(order >0)
    {
        for(int i=0; i<_sub_clst_pos.size(); ++i)
        {
            CluserPostion p ( order-1, _sub_clst_pos[i] );
            res.insert( p );
        }
    }
    return res;
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
    size_t sz1=idx1._index.size(); size_t sz2=idx2._index.size();
    if( sz1 == sz2 )
        for(int i=0; i<sz1; ++i)
        { if(idx1._index[i] != idx2._index[i]) return idx1._index[i] < idx2._index[i]; }
    else
        return (sz1 < sz2);
    return 0;// idx1 equals to idx2
}

ostream&  operator << (ostream& outs, const cClusterIndex& idx)
{
    outs << "[" ;
    for(int i=0; i<idx._index.size(); ++i)
    {
        outs << idx._index(i);
        if(i<idx._index.size()-1)
            outs <<", ";
    }
    outs << "];\t";
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
{ //LOG(INFO) << "Constructor of cSpinGrouping with connextion_matrix.";
    _connection_matrix=connection_matrix;
}

cSpinGrouping::~cSpinGrouping()
{ //LOG(INFO) << "Default destructor of cSpinGrouping";
}

sp_mat cSpinGrouping::index2subgraph(int order)
{
    size_t nClst=_cluster_index_list[order].size();
    mat res=zeros(nClst, _nspin);

//    int i=0;
//    for(cClusterIndex vIdx:_cluster_index_list[order])
//    {
//        res.row(i) = vIdx.get_array(_nspin);
//        i++;
//    }
    int i=0;
    FIX_ORDER_INDEX_SET clst_set = _cluster_index_list[order];
    for(set<cClusterIndex>::iterator pos=clst_set.begin(); pos!=clst_set.end(); ++pos)
    {
        cClusterIndex vIdx = *pos;
        res.row(i) = vIdx.get_array(_nspin);
        i++;
    }
    return conv_to<sp_mat>::from(res);
}

void cSpinGrouping::subgraph2index(const sp_mat& subgraph, const vector<int> sub_pos_list)
{
    for(int i=0; i<subgraph.n_rows; ++i)
    {
        //cout << "\r" << setw(6) <<  i+1 << "/" << subgraph.n_rows 
             //<< " subgraphs are inserted.";
        mat r(subgraph.row(i));  uvec nz_r = find(r);  size_t order = nz_r.size()-1;
        cClusterIndex cIdx( nz_r );
        pair<FIX_ORDER_INDEX_SET::iterator, bool> pos = _cluster_index_list[ order ].insert(cIdx);
        if( order > 0)
            pos.first->appendSubClstPos( sub_pos_list[i] );
    }
    //cout <<endl;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinDepthFirstPathTracing
cDepthFirstPathTracing::cDepthFirstPathTracing()
{ //LOG(INFO) << "Default constructor: cDepthFirstPathTracing.";
}

cDepthFirstPathTracing::cDepthFirstPathTracing(const sp_mat&  connection_matrix, size_t maxOrder)
{
    _max_order = maxOrder;
    _nspin     = connection_matrix.n_cols;
    _connection_matrix=connection_matrix;
    vector<int> empty (0);
    subgraph2index( speye(_nspin, _nspin), empty );
}

cDepthFirstPathTracing::cDepthFirstPathTracing(const sp_mat&  connection_matrix, size_t maxOrder, const mat& init)
{
    _max_order = maxOrder;
    _nspin     = connection_matrix.n_cols;
    _connection_matrix=connection_matrix;
    vector<int> empty (0);
    sp_mat init_spmat=conv_to<sp_mat>::from( init );
    subgraph2index( init_spmat, empty );
}

cDepthFirstPathTracing::~cDepthFirstPathTracing()
{ //LOG(INFO) << "Default destructor: cDepthFirstPathTracing.";
}

void cDepthFirstPathTracing::generate()
{
    sp_mat subgraph = index2subgraph(0);

    for( int i = 1; i < _max_order; ++i)
    {
        sp_mat neighbor = subgraph*_connection_matrix;
        pair<sp_mat, vector<int> > growth_res = subgraph_growth(subgraph, neighbor, i);
        sp_mat new_subgraph = growth_res.first;
        vector<int> sub_pos = growth_res.second;
        subgraph2index(new_subgraph, sub_pos);
        subgraph= index2subgraph(i);
    }
}

pair<sp_mat, vector<int> >  cDepthFirstPathTracing::subgraph_growth(const sp_mat& subgraph, const sp_mat& neighbor, int subgraph_order)
//sp_mat cDepthFirstPathTracing::subgraph_growth(const sp_mat& subgraph, const sp_mat& neighbor, int subgraph_order)
{
    mat new_subgraph=zeros( sum(nonzeros(neighbor)), _nspin);
    vector<int> sub_pos_list;

    int nGen=0;
    for(int i=0; i<subgraph.n_rows; ++i)
    {
        //cout << "\r"<< i+1 <<  "/" << subgraph.n_rows
             //<< " clusters are generated of spin order= " << subgraph_order+1 << "\t";

        mat parent_row(subgraph.row(i));
        mat r(neighbor.row(i));    uvec candidate = find(r);

        //size_t row_idx=0;
        for(int j=0; j<candidate.size(); ++j)
        {
            if(parent_row(candidate(j))==0)
            {
                new_subgraph.row(nGen)=parent_row;
                new_subgraph(nGen, candidate(j))=1;
                sub_pos_list.push_back( i );
                nGen++;
            }
        }
        //cout.flush();
    }
    //cout << endl;

    sp_mat res_mat=conv_to<sp_mat>::from( new_subgraph.rows(0, nGen-1) );

    pair<sp_mat, vector<int> > res (res_mat, sub_pos_list);
    //return res_mat;
    return res;
}

//}}}
////////////////////////////////////////////////////////////////////////////////
