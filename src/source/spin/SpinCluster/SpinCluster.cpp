#include "include/spin/SpinCluster.h"
//#include "include/spin/SpinClusterAlgorithm.h"

////////////////////////////////////////////////////////////////////
//{{{ cClusterIndex
cClusterIndex::cClusterIndex()
{ LOG(INFO) << "Default constructor of cClusterIndex.";
}

cClusterIndex::cClusterIndex(const uvec& idx)
{
    _index = idx;
    sort(_index.begin(), _index.end()); }

cClusterIndex::~cClusterIndex()
{ LOG(INFO) << "Default destructor of cClusterIndex.";
}

mat cClusterIndex::get_array(size_t nspin)
{
    mat idx_array=zeros(1, nspin);
    int nnz = _index.size();

    for(int i=0; i<nnz; ++i)
        idx_array[_index[i]]=1;
    return idx_array;
}

bool operator == (const cClusterIndex& idx1, const cClusterIndex& idx2)
{
    if(idx1._index.size() == idx2._index.size())
        for(int i=0; i<idx1._index.size(); ++i)
        { if(idx1._index[i] != idx2._index[i]) return 0; }
    else
        return 0;
    return 1;
}


bool operator < (const cClusterIndex& idx1, const cClusterIndex& idx2)
{
    int sz1=idx1._index.size(); int sz2=idx2._index.size();
    if( sz1 == sz2 )
        for(int i=0; i<sz1; ++i)
        { if(idx1._index[i] != idx2._index[i]) return idx1._index[i] < idx2._index[i]; }
    else
        return (sz1 < sz2);
    return 0;// idx1 equals to idx2
}

ostream&  operator << (ostream& outs, const cClusterIndex& idx)
{
    for(auto it=idx._index.begin(); it !=idx._index.end(); ++it)
    {
        outs << *it;
        if(next(it) != idx._index.end())
            outs << ", ";
    }
    return outs;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////
//{{{ cSpinCluster
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
/// This function calls the 'generate' method of the grouping algorithm.
    _grouping->generate();
    _cluster_index_list = _grouping->get_cluster_index();
}

ostream&  operator << (ostream& outs, const cSpinCluster& clst)
{
/// Operator << is reloaded to display the cluster index one by one.
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
//}}}
////////////////////////////////////////////////////////////////////
