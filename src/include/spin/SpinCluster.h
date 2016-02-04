#ifndef SPINCLUSTER_H
#define SPINCLUSTER_H

#include <vector>
#include <set>
#include "include/easylogging++.h"
#include "include/spin/Spin.h"
#include "include/spin/SpinCollection.h"
#include "include/spin/SpinClusterAlgorithm.h"

/// \addtogroup SpinList
/// @{

/// \defgroup SpinCluster SpinCluster
/// @{
typedef vector<umat> clusterTable; // nWorder <umat>

struct MPI_Cluster_Data
{
    int nWorker;
    int nOrder;
    umat jobTable; // nOrder * nWorker 

    vector<int>  clusterNumList; // nOrder <int> 
    vector<clusterTable> clusterData; // nOrder <clusterTable>
};

////////////////////////////////////////////////////////////////////
//{{{ cSpinCluster
/// This class generates spin clusters with a given grouping algorithm.
///
class cSpinCluster
{
public:
    cSpinCluster();
    cSpinCluster(const cSpinCollection& sc, cSpinGrouping * grouping);
    cSpinCluster(const cSpinCollection& sc, const uvec& clstLength, const vector<umat>& clstMatList);
    ~cSpinCluster();

    void make();
    void MPI_partition(int nWorker);

    CLST_IDX_LIST getClusterIndex() const {return _cluster_index_list;};
    umat          getClusterIndex(size_t order) const ;
    cClusterIndex getClusterIndex(const CluserPostion& pos) const {return getClusterIndex(pos.first, pos.second);};
    cClusterIndex getClusterIndex(size_t order, size_t index) const ;
    vector<cSPIN> getCluster(size_t order, size_t index) const ;
    size_t        getMaxOrder() const {return _max_order;};
    size_t        getClusterNum(int order) const {return _cluster_index_list[order].size();};
    set<CluserPostion > getSubClusters(size_t order, size_t index) const;

    uvec          getMPI_ClusterLength(int worker_id) const {return _data.jobTable.col(worker_id);};
    vector<umat>  getMPI_Cluster(int worker_id);
    CluserPostion getMPI_ClusterSize(int cce_order, int worker_id) const;

    friend ostream&  operator << (ostream& outs, const cSpinCluster& clst);
private:
    size_t           _max_order;
    cSpinGrouping *  _grouping;
    CLST_IDX_LIST    _cluster_index_list;
    cSpinCollection  _spin_collection;
    MPI_Cluster_Data _data;
};
//}}}
////////////////////////////////////////////////////////////////////

/// @}
/// @}
#endif
