#include "include/spin/SpinCluster.h"
#include "include/spin/SpinClusterAlgorithm.h"
////////////////////////////////////////////////////////////////////////////////
//{{{ cUniformBathOnLattice
class cUniformBathOnLattice:public cSpinGrouping
{
public:
    cUniformBathOnLattice(){};
    cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins);
    ~cUniformBathOnLattice(){};

    void generate();
private:
    void generate_primitive_clusters();

    int                    _bath_center_index;
    cSpinCollection        _bath_spins;
    cDepthFirstPathTracing _primitive_dfpt;
    cSpinCluster           _primitive_spin_clusters;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
