#ifndef SPINCLUSTER_H
#define SPINCLUSTER_H

#include <vector>
#include <set>
#include "include/spin/Spin.h"
#include "include/spin/SpinGrouping.h"

class cSpinCluster
{
public:
    cSpinCluster(cSpinGrouping * grouping);
    ~cSpinCluster();

    void make();
    CLST getClusterIndex(){return _cluster_index_list;};
private:
    cSpinGrouping * _grouping;
    CLST _cluster_index_list;
};
#endif
