#ifndef SPINCLUSTERFROMLATTICE_H
#define SPINCLUSTERFROMLATTICE_H

#include "include/spin/SpinCluster.h"
#include "include/spin/SpinClusterAlgorithm.h"

typedef vector<size_t> SubPosLst;
typedef vector< SubPosLst > SubPosLst_FixOrder;
////////////////////////////////////////////////////////////////////////////////
//{{{ cUniformBathOnLattice
class cUniformBathOnLattice:public cSpinGrouping
{
public:
    cUniformBathOnLattice(){};
    cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins, const Lattice& lattice, int root_range_idx);
    cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins, const Lattice& lattice, const imat& root_range);
    ~cUniformBathOnLattice(){};

    void generate();
    int  getPrimitiveClusterNumber(int corner, int order) const { return _primitive_spin_clusters[corner].getClusterNum(order);}


private:
    Lattice                _lattice;
    Lattice                _root_lattice;
    vector<int>            _center;
    int                    _unit_cell_num;
    int                    _atom_num_in_cell;
    imat                   _root_range;
    vector<int>            _root_center;
    vector<int>            _root_index_list;

    cSpinCollection        _bath_spins;
    vector<cSPIN>          _spin_list;

    vector<cSpinCluster>   _primitive_spin_clusters;
    umat                   _primitive_cluster_size;
    umat                   _primitive_cumsum_size;
    vector<umat>           _primitive_cluster_mat;

    vector<int>                 _total_cluster_number;
    vector< SubPosLst_FixOrder> _sub_pos;


    void      generate_primitive_clusters();
    void      generate_sub_primitive_position();
    void      generate_cluster_index_list();
    int       locate_primitive_sub_clusters(const urowvec& v);
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
