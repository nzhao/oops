#include "include/spin/SpinCollection.h"
#include "include/spin/SpinCluster.h"
//#include "include/spin/SpinClusterAlgorithm.h"



////////////////////////////////////////////////////////////////////
//{{{ cSpinCluster
cSpinCluster::cSpinCluster()
{ LOG(INFO) << "Defaul constructor: cSpinCluster.";
}
cSpinCluster::cSpinCluster(const cSpinCollection& sc, cSpinGrouping * grouping)
{
    _grouping = grouping;
    _spin_collection = sc;
}

cSpinCluster::~cSpinCluster()
{
    if (!_grouping) delete _grouping;
}

void cSpinCluster::make()
{
/// This function calls the 'generate' method of the grouping algorithm.
    _grouping->generate();
    _cluster_index_list = _grouping->get_cluster_index();
}

cClusterIndex cSpinCluster::getClusterIndex(int order, int index) const
{
    FIX_ORDER_INDEX_SET::iterator it = _cluster_index_list[order].begin();
    advance(it, index);
    return *it;
}

vector<cSPIN> cSpinCluster::getCluster(int order, int index) const
{
    cClusterIndex clst = getClusterIndex(order, index);
    return _spin_collection.getSpinList(clst);
}

ostream&  operator << (ostream& outs, const cSpinCluster& clst)
{
/// Operator << is reloaded to display the cluster index one by one.
    int i, j, tot; i=1; j=1; tot=0;
    
    outs << "Total Order = " << clst._cluster_index_list.size() << endl;
//    for(auto clst_set: clst._cluster_index_list)
//    {
//        j=1; tot += clst_set.size();
//        if(clst_set.size() > 0)
//        {
//            outs << "Cluster Order = " << i << ": Number = " << clst_set.size() << ": " << endl;
//            for(auto idx: clst_set)
//            {
//                outs << j << ": " <<  idx << endl;
//                j++;
//            }
//            outs << endl;
//        }
//        i++;
//    }
    for(int order=0; order<clst._cluster_index_list.size(); ++order)
    {
        FIX_ORDER_INDEX_SET clst_set = clst._cluster_index_list[order];

        j=1; tot += clst_set.size();
        if(clst_set.size() > 0)
        {
            outs << "Cluster Order = " << i << ": Number = " << clst_set.size() << ": " << endl;
            for(set<cClusterIndex>::iterator pos=clst_set.begin(); pos!=clst_set.end(); ++pos)
            {
                cClusterIndex vIdx = *pos;
                outs << j << ": " <<  vIdx << endl;
                j++;
            }
            outs << endl;
        }
        i++;
    }
    outs << tot << " clusters are generated." << endl;
    return outs;
}
//}}}
////////////////////////////////////////////////////////////////////
