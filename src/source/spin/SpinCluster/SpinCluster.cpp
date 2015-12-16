#include "include/spin/SpinCluster.h"
#include "include/spin/SpinClusterAlgorithm.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
// cSpinCluster
cSpinCluster::cSpinCluster(cSpinGrouping * grouping)
{
    _grouping = grouping;
}

cSpinCluster::~cSpinCluster()
{
    if (!_grouping) delete _grouping;
}

void cSpinCluster::make()
{
    _grouping->generate();
    _cluster_index_list = _grouping->get_cluster_index();
}

ostream&  operator << (ostream& outs, const cSpinCluster& clst)
{
    for(auto clst_set: clst._cluster_index_list)
        for(auto idx: clst_set)
            outs << idx << endl;
    return outs;
}
