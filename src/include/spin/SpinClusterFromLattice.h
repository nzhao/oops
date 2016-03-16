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
    cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins, const Lattice& lattice);
    ~cUniformBathOnLattice(){};

    void        generate();
    int         getPrimitiveClusterNumber(int corner, int order) const { return _primitive_spin_clusters[corner].getClusterNum(order);}
    void        generate_cluster_index_list();
    //{{{
    //SubPrPosLst find_sub_position(const primitive_position& pos);

    //uvec        getPrimitiveClusterIndex(const primitive_position& pos) const;
    //uvec        getPrimitiveClusterIndex(int corner, int order, int idx) const {return _primitive_spin_clusters[corner].getClusterIndex(order, idx).getIndex();}
    //vector<vec> getPrimitiveClusterCoord(const primitive_position& pos) const { return _primitive_spin_clusters[pos.getCorner()].getClusterCoord( pos.getOrder(), pos.getIndex() );}

    //int         getGlobalPosition(int unit_cell_index, const primitive_position& pos);
    //int         getGlobalClusterNumber(int order) const {return _lattice.getUnitCellNumber() * _primitive_cumsum_size(_atom_num_in_cell, order);}
    //uvec        getClusterIndex(int unit_cell, int corner_idx, int order, int idx);
    //}}}


private:
    Lattice                _lattice;
    vector<int>            _center;
    int                    _unit_cell_num;
    int                    _atom_num_in_cell;

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
    
    int       locate_primitive_sub_clusters(const urowvec& v);

    //CLST_IDX_LIST _cluster_index_list;
    //{{{
    //pair<int, int>  match_pos(const urowvec& idx);
    //primitive_position compute_sub_pos(const primitive_position& pos, size_t remove_k);
    //vector<vec> relative_coordinates(const primitive_position& pos);
    //pair< vector<int>, vector<vec> > remove_k(const primitive_position& pos, size_t remove_k);
    //double patten_diff(const vector<vec>& v1, const vector<vec>& v2);
    //}}}
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
