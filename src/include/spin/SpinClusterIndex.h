#ifndef SPINCLUSTERINDEX_H
#define SPINCLUSTERINDEX_H

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <armadillo>
#include "include/spin/Spin.h"

typedef pair<size_t, size_t> ClusterPostion;

////////////////////////////////////////////////////////////////////////////////
//{{{ cClusterIndex
/// This class defines index list of a cluster.
///
class cClusterIndex
{
public:
    cClusterIndex();
    cClusterIndex(const uvec& idx);
    ~cClusterIndex();

    mat get_array(size_t nspin);
    uvec getIndex() const {return _index;};
    size_t getNum() const {return _spin_num;};
    size_t getOrder() const {return _spin_num-1;};
    set< ClusterPostion > getSubClstPos() const;

    void appendSubClstPos(int pos) const  {_sub_clst_pos.push_back( pos );};
    void setSubClstPos(vector<size_t> pos_lst) const {_sub_clst_pos = pos_lst;};

    friend bool operator == (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend bool operator < (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend ostream&  operator << (ostream& outs, const cClusterIndex& idx);
private:
    uvec _index;
    size_t _spin_num;
    mutable vector<size_t> _sub_clst_pos;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
