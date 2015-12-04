#ifndef SPINCLUSTER_H
#define SPINCLUSTER_H

#include <vector>
#include "include/spin/Spin.h"
#include "include/spin/SpinGathering.h"

typedef vector<vector<int>> CLST_IDX;

class cSpinCluster
{
public:
    cSpinCluster(cSpinGathering * gathering);
    ~cSpinCluster();

    void make();
    CLST_IDX getClusterIndex(){return _groupling->get_cluster_index();};
private:
    cSpinGathering * _gathering;
};
#endif
