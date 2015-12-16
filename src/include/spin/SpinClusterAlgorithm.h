#ifndef SPINGROUPING_H
#define SPINGROUPING_H
#define MAX_CLUSTER_ORDER  10

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <armadillo>
#include "include/spin/Spin.h"

using namespace std;
using namespace arma;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cClusterIndex
class cClusterIndex
{
public:
    cClusterIndex();
    cClusterIndex(const uvec& idx);
    ~cClusterIndex();

    sp_mat get_array(size_t nspin);

    friend bool operator == (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend bool operator < (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend ostream&  operator << (ostream& outs, const cClusterIndex& idx);
private:
    uvec _index;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinGrouping

typedef set<cClusterIndex>          FIX_ORDER_INDEX_SET;
typedef vector<FIX_ORDER_INDEX_SET> CLST_IDX_LIST;

class cSpinGrouping
{
public:
    cSpinGrouping(const sp_mat& connection_matrix);
    cSpinGrouping();
    virtual ~cSpinGrouping();
    virtual void generate()=0;

    CLST_IDX_LIST& get_cluster_index() {return _cluster_index_list;};
protected:
    size_t        _nspin;
    sp_mat        _connection_matrix;
    CLST_IDX_LIST _cluster_index_list;

    void subgraph2index(const sp_mat& subgraph);
    sp_mat index2subgraph(int order);
private:
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSourceFromFile
class cDepthFirstPathTracing:public cSpinGrouping
{
public:
    cDepthFirstPathTracing();
    cDepthFirstPathTracing(const sp_mat& connection_matrix, size_t maxOrder);
    virtual ~cDepthFirstPathTracing();

    void generate();

private:
    size_t _max_order;
    sp_mat subgraph_growth(const sp_mat& subgraph, const sp_mat& neighbor, int subgraph_order);

};
#endif
