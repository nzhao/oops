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
    CLST_IDX_LIST getClusterIndex(){return _cluster_index_list;};
    umat          getClusterIndex(int order) const ;
    cClusterIndex getClusterIndex(int order, int index) const ;
    vector<cSPIN> getCluster(int order, int index) const ;

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
