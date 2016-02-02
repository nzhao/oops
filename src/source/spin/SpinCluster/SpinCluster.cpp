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

cClusterIndex cSpinCluster::getClusterIndex(size_t order, size_t index) const
{
    FIX_ORDER_INDEX_SET::iterator it = _cluster_index_list[order].begin();
    advance(it, index);
    return *it;
}

umat cSpinCluster::getClusterIndex(size_t order) const
{
    umat res = zeros<umat> (getClusterNum(order), order+1);

    FIX_ORDER_INDEX_SET::iterator it;
    const FIX_ORDER_INDEX_SET& clusters = _cluster_index_list[order];
    for(it = clusters.begin(); it != clusters.end(); ++it)
        res.row( distance(clusters.begin(), it ) ) = trans( it->getIndex() );
    return res;
}

MPI_Cluster_Data cSpinCluster::MPI_partition(int nWorker)
{
    MPI_Cluster_Data res;
    res.nWorker = nWorker;
    res.nOrder = getMaxOrder();

    res.jobTable = umat(res.nOrder+1, res.nWorker, fill::zeros);
    for(int order_i = 0; order_i<res.nOrder; ++order_i)
    {
        int clstNum = getClusterNum(order_i);
        res.clusterNumList.push_back( clstNum );

        umat full_clst_idx = getClusterIndex( order_i );
        clusterTable clst_tb_i;

        int row1 = 0; int row2 = 0;
        int q = clstNum/nWorker; int r = clstNum % nWorker;
        for(int wk_id = 0; wk_id<nWorker; ++wk_id)
        {
            int jobs = wk_id < r ? q + 1: q;
            res.jobTable(order_i, wk_id) = jobs;

            row2 = row1 + jobs;
            clst_tb_i.push_back( full_clst_idx.rows(row1, row2) );
            row1 = row2;
        }
        res.clusterData.push_back(clst_tb_i);
    }
    return res;
}

vector<cSPIN> cSpinCluster::getCluster(size_t order, size_t index) const
{
    cClusterIndex clst = getClusterIndex(order, index);
    return _spin_collection.getSpinList(clst);
}

ostream&  operator << (ostream& outs, const cSpinCluster& clst)
{/*{{{*/
/// Operator << is reloaded to display the cluster index one by one.
    int i, j, tot; i=1; j=1; tot=0;
    
    outs << "Total Order = " << clst._cluster_index_list.size() << endl;
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
}/*}}}*/
//}}}
////////////////////////////////////////////////////////////////////
