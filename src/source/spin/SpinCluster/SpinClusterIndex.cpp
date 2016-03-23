#include "include/spin/SpinClusterIndex.h"

////////////////////////////////////////////////////////////////////
//{{{ cClusterIndex
cClusterIndex::cClusterIndex()
{ //LOG(INFO) << "Default constructor of cClusterIndex.";
}

cClusterIndex::cClusterIndex(const uvec& idx)
{
    _index = idx;
    sort(_index.begin(), _index.end());
    _spin_num = idx.n_elem;
}

cClusterIndex::~cClusterIndex()
{ //LOG(INFO) << "Default destructor of cClusterIndex.";
}

mat cClusterIndex::get_array(size_t nspin)
{
    mat idx_array=zeros(1, nspin);
    size_t nnz = _index.size();

    for(int i=0; i<nnz; ++i)
        idx_array[_index[i]]=1;
    return idx_array;
}
set< ClusterPostion > cClusterIndex::getSubClstPos() const
{
    set< ClusterPostion > res;
    size_t order = getOrder();
    if(order >0)
    {
        for(int i=0; i<_sub_clst_pos.size(); ++i)
        {
            ClusterPostion p ( order-1, _sub_clst_pos[i] );
            res.insert( p );
        }
    }
    return res;
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
    size_t sz1=idx1._index.size(); size_t sz2=idx2._index.size();
    if( sz1 == sz2 )
        for(int i=0; i<sz1; ++i)
        { if(idx1._index[i] != idx2._index[i]) return idx1._index[i] < idx2._index[i]; }
    else
        return (sz1 < sz2);
    return 0;// idx1 equals to idx2
}

ostream&  operator << (ostream& outs, const cClusterIndex& idx)
{
    outs << "[" ;
    for(int i=0; i<idx._index.size(); ++i)
    {
        outs << idx._index(i);
        if(i<idx._index.size()-1)
            outs <<", ";
    }
    outs << "];";
    //for(int i=0; i<idx._sub_clst_pos.size(); ++i)
    //{
        //outs << idx._sub_clst_pos[i];
        //if(i<idx._sub_clst_pos.size()-1)
            //outs <<", ";
    //}
    //outs << "};\t";
    return outs;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
