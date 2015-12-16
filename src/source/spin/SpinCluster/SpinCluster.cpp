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
    int i, j, tot; i=1; j=1; tot=0;
    
    outs << "Total Order = " << clst._cluster_index_list.size() << endl;
    for(auto clst_set: clst._cluster_index_list)
    {
        j=1; tot += clst_set.size();
        if(clst_set.size() > 0)
        {
            outs << "Cluster Order = " << i << ": Number = " << clst_set.size() << ": " << endl;
            for(auto idx: clst_set)
            {
                outs << j << ": " <<  idx << endl;
                j++;
            }
            outs << endl;
        }
        i++;
    }
    outs << tot << " clusters are generated." << endl;
    return outs;
}
