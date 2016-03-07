#ifndef SPINGROUPING_H
#define SPINGROUPING_H
#define MAX_CLUSTER_ORDER  10

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <armadillo>
#include "include/easylogging++.h"
#include "include/spin/Spin.h"
//#include "include/spin/SpinCluster.h"

using namespace std;
using namespace arma;

typedef pair<size_t, size_t> CluserPostion;

 
/// \addtogroup SpinList
/// @{

/// \addtogroup SpinCluster
/// @{

////////////////////////////////////////////////////////////////////////////////
//{{{ cClusterIndex
/// This class defines index list of a cluster.
///
class cClusterIndex
{
public:
    cClusterIndex();
    cClusterIndex(const uvec& idx);
    ~cClusterIndex();

    mat get_array(size_t nspin);
    uvec getIndex() const {return _index;};
    size_t getNum() const {return _spin_num;};
    size_t getOrder() const {return _spin_num-1;};
    set< CluserPostion > getSubClstPos() const;

    void appendSubClstPos(int pos) const  {_sub_clst_pos.push_back( pos );};

    friend bool operator == (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend bool operator < (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend ostream&  operator << (ostream& outs, const cClusterIndex& idx);
private:
    uvec _index;
    size_t _spin_num;
    mutable vector<size_t> _sub_clst_pos;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// \defgroup SpinGrouping
/// @{

typedef set<cClusterIndex>          FIX_ORDER_INDEX_SET;
typedef vector<FIX_ORDER_INDEX_SET> CLST_IDX_LIST;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinGrouping

/// This is an abstract class of spin grouping algorithm.
///
class cSpinGrouping
{
public:
    cSpinGrouping(const sp_mat& connection_matrix);
    cSpinGrouping();
    virtual ~cSpinGrouping();
    virtual void generate()=0;

    size_t         getMaxOrder() const {return _max_order;};
    CLST_IDX_LIST& get_cluster_index() {return _cluster_index_list;};
protected:
    size_t        _nspin;
    size_t        _max_order;
    sp_mat        _connection_matrix;
    CLST_IDX_LIST _cluster_index_list;

    void subgraph2index(const sp_mat& subgraph, const vector<int> sub_pos_list);
    sp_mat index2subgraph(int order);
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ cDepthFirstPathTracing
/// This class implements a concreate grouping algorithm by 'depth-first search'. Ref. Spinach. 
///
class cDepthFirstPathTracing:public cSpinGrouping
{
public:
    cDepthFirstPathTracing();
    cDepthFirstPathTracing(const sp_mat& connection_matrix, size_t maxOrder);
    cDepthFirstPathTracing(const sp_mat&  connection_matrix, size_t maxOrder, const sp_mat& init_graph);
    virtual ~cDepthFirstPathTracing();

    void generate();

private:
    pair<sp_mat, vector<int> >  subgraph_growth(const sp_mat& subgraph, const sp_mat& neighbor, int subgraph_order);
    //sp_mat remove_repeat(sp_mat subgraph, int subgraph_order);

};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}
/// @}
/// @}
#endif
