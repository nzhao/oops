#ifndef SPINCLUSTER_H
#define SPINCLUSTER_H

#include <vector>
#include <set>
#include "include/easylogging++.h"
#include "include/spin/Spin.h"
#include "include/spin/SpinClusterAlgorithm.h"

/// \addtogroup SpinList
/// @{

/// \defgroup SpinCluster SpinCluster
/// @{
typedef vector<umat> clusterTable; // nWorder <uma>

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
    ~cSpinCluster();

    void make();
    MPI_Cluster_Data MPI_partition(int nWorker);

    CLST_IDX_LIST getClusterIndex() const {return _cluster_index_list;};
    umat          getClusterIndex(size_t order) const ;
    cClusterIndex getClusterIndex(size_t order, size_t index) const ;
    vector<cSPIN> getCluster(size_t order, size_t index) const ;
    size_t        getMaxOrder() const {return _grouping->getMaxOrder();};
    size_t        getClusterNum(int order) const {return _cluster_index_list[order].size();};

    friend ostream&  operator << (ostream& outs, const cSpinCluster& clst);
private:
    cSpinGrouping * _grouping;
    CLST_IDX_LIST _cluster_index_list;
    cSpinCollection _spin_collection;
};
//}}}
////////////////////////////////////////////////////////////////////

/// @}
/// @}
#endif
