#include "include/spin/SpinCollection.h"
#include "include/spin/SpinCluster.h"
//#include "include/spin/SpinClusterAlgorithm.h"



////////////////////////////////////////////////////////////////////
//{{{ cSpinCluster
cSpinCluster::cSpinCluster()
{ //LOG(INFO) << "Defaul constructor: cSpinCluster.";
}
cSpinCluster::cSpinCluster(const cSpinCollection& sc, cSpinGrouping * grouping)
{
    _grouping = grouping;
    _spin_collection = sc;
    _max_order = grouping->getMaxOrder();
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
    _cluster_index_mat = _grouping->get_cluster_index_mat();
    _sub_cluster_position_valid = true;
}

cSpinCluster::cSpinCluster(const cSpinCluster& clst)
{
    _spin_collection = clst._spin_collection;
    _max_order = clst._max_order;//clst.getMaxOrder();
    _cluster_index_list = clst._cluster_index_list;//clst.getClusterIndex();
}

cSpinCluster::cSpinCluster(const cSpinCollection& sc, const uvec& clstLength, const vector<umat>& clstMatList)
{
    _spin_collection = sc;
    _max_order = clstMatList.size();
    for(int i=0; i<clstMatList.size(); ++i)
    {
        umat fix_order_mat = clstMatList[i];
        FIX_ORDER_INDEX_SET fix_order_set;
        for(int j=0; j<fix_order_mat.n_rows; ++j)
        {
            uvec v = trans( fix_order_mat.row(j) );
            fix_order_set.insert( cClusterIndex(v) );
        }
        _cluster_index_list.push_back(fix_order_set);
    }
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

    if( order < _cluster_index_mat.size() )
    {
        cout << _cluster_index_mat.size() << endl;
        res = _cluster_index_mat[order];
    }
    else
    {
        FIX_ORDER_INDEX_SET::iterator it;
        const FIX_ORDER_INDEX_SET& clusters = _cluster_index_list[order];
        for(it = clusters.begin(); it != clusters.end(); ++it)
            res.row( distance(clusters.begin(), it ) ) = trans( it->getIndex() );
    }
    return res;
}
set<ClusterPostion > cSpinCluster::getSubClusters(size_t order, size_t index) const
{
    //cout << "{ " << order << ", " << index << " }." << endl;
    cClusterIndex clst = getClusterIndex(order, index);
    //cout << clst << endl;
    set<ClusterPostion > sub_pos = clst.getSubClstPos();
    set<ClusterPostion >::iterator it;
    for(it=sub_pos.begin(); it!=sub_pos.end(); ++it)
    {
        //cout << "\t  it = {" << it->first << ", " << it->second << "} " << endl;
        set<ClusterPostion > sub_sub_pos = getSubClusters(it->first, it->second);//sub_cluster.getSubClstPos();
        set<ClusterPostion >::iterator sub_it;
        for(sub_it=sub_sub_pos.begin(); sub_it!=sub_sub_pos.end(); ++sub_it)
            sub_pos.insert( *sub_it );
    }
    return sub_pos;
}

void cSpinCluster::MPI_partition(int nWorker)
{
    _data.nWorker = nWorker;
    _data.nOrder = getMaxOrder();

    _data.jobTable = umat(_data.nOrder, _data.nWorker, fill::zeros);
    for(int order_i = 0; order_i<_data.nOrder; ++order_i)
    {
        cout << "partitioning order = " << order_i << endl;
        int clstNum = getClusterNum(order_i);
        _data.clusterNumList.push_back( clstNum );

        umat full_clst_idx = getClusterIndex( order_i );
        //cout << full_clst_idx << endl;
        clusterTable clst_tb_i;

        int row1 = 0; int row2 = 0;
        int q = clstNum/nWorker; int r = clstNum % nWorker;

        for(int wk_id = 0; wk_id<nWorker; ++wk_id)
        {
            cout << "worker_id = " << wk_id << endl;
            int jobs = wk_id < r ? q + 1: q;
            _data.jobTable(order_i, wk_id) = jobs;

            row2 = row1 + jobs;
            
            clst_tb_i.push_back( full_clst_idx.rows(row1, row2-1) );
            row1 = row2;
        }
        _data.clusterData.push_back(clst_tb_i);
    }
}

vector<umat> cSpinCluster::getMPI_Cluster(int worker_id)
{
    vector<umat> res;
    for(int i=0; i<getMaxOrder(); ++i)
        res.push_back( _data.clusterData[i][worker_id] );
    return res;
}

ClusterPostion cSpinCluster::getMPI_ClusterSize(int cce_order, int worker_id) const
{
    size_t pos1, pos2;
    pos1 = 0;
    for(int id=0; id<worker_id; ++id)
        pos1 += _data.jobTable(cce_order, id); 
    pos2 = pos1 + _data.jobTable(cce_order, worker_id);

    ClusterPostion res(pos1, pos2);
    return res;
}

vector<cSPIN> cSpinCluster::getCluster(size_t order, size_t index) const
{
    cClusterIndex clst = getClusterIndex(order, index);
    return _spin_collection.getSpinList(clst);
}

vector<vec> cSpinCluster::getClusterCoord(size_t order, size_t index) const
{
    vector<vec> coord_list;
    vector<cSPIN> sl = getCluster(order, index);
    for(int i=0; i<sl.size(); ++i)
        coord_list.push_back( sl[i].get_coordinate() );
    return coord_list;
}

ostream&  operator << (ostream& outs, const cSpinCluster& clst)
{/*{{{*/
/// Operator << is reloaded to display the cluster index one by one.
    int i, j, tot; i=0; j=0; tot=0;
    
    outs << "Total Order = " << clst._cluster_index_list.size() << endl;
    for(int order=0; order<clst._cluster_index_list.size(); ++order)
    {
        FIX_ORDER_INDEX_SET clst_set = clst._cluster_index_list[order];

        j=0; tot += clst_set.size();
        if(clst_set.size() > 0)
        {
            outs << "Cluster Order = " << i << ": Number = " << clst_set.size() << ": " << endl;
            if(clst_set.size() > 100)
                cout << "first 100 clusters are displayed." << endl;
            for(set<cClusterIndex>::iterator pos=clst_set.begin(); pos!=clst_set.end() && j<100; ++pos)
            {
                cClusterIndex vIdx = *pos;
                outs << "{ " << order << ", " << j << " } = "  <<  vIdx << "\t" ;

                if(clst._sub_cluster_position_valid)
                {
                    set<ClusterPostion > sub_pos = clst.getSubClusters(order, j);//vIdx.getSubClstPos();
                    set<ClusterPostion >::iterator it;
                    for(it=sub_pos.begin(); it!=sub_pos.end(); ++it)
                        outs << "{ " << it->first << ", " << it->second << " }\t" ;
                }
                outs << endl;
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
