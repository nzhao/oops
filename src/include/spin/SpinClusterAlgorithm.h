#ifndef SPINGROUPING_H
#define SPINGROUPING_H
#define MAX_CLUSTER_ORDER  10

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <armadillo>
#include "include/spin/Spin.h"
//#include "include/misc/lattice.h"
#include "include/spin/SpinCollection.h"
#include "include/spin/SpinClusterIndex.h"

using namespace std;
using namespace arma;

 
/// \addtogroup SpinList
/// @{

/// \addtogroup SpinCluster
/// @{


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
    vector<umat>   get_cluster_index_mat() const {return _cluster_index_mat;};

protected:
    size_t        _nspin;
    size_t        _max_order;
    sp_mat        _connection_matrix;
    CLST_IDX_LIST _cluster_index_list;
    vector<umat>  _cluster_index_mat;

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
    cDepthFirstPathTracing(const sp_mat& connection_matrix, size_t maxOrder, const mat& init);
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
