#include "include/spin/SpinCluster.h"
#include "include/spin/SpinClusterAlgorithm.h"
class  primitive_position
{
public:
    primitive_position(){};
    primitive_position(int c, int o, int i): _corner(c), _order(o), _index(i){};
    ~primitive_position(){};

    int getCorner() const {return _corner;}
    int getOrder() const {return _order;}
    int getIndex() const {return _index;}
    friend ostream& operator << (ostream& outs, const primitive_position& pos);
private:
    int _corner;
    int _order;
    int _index;
};

////////////////////////////////////////////////////////////////////////////////
//{{{ cUniformBathOnLattice
class cUniformBathOnLattice:public cSpinGrouping
{
public:
    cUniformBathOnLattice(){};
    cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins, const Lattice& lattice);
    ~cUniformBathOnLattice(){};

    void                       generate();
    void                       generate_cluster_index_list();
    vector<primitive_position> find_sub_position(const primitive_position& pos);

    uvec        getPrimitiveClusterIndex(const primitive_position& pos) const;
    uvec        getPrimitiveClusterIndex(int corner, int order, int idx) const {return _primitive_spin_clusters[corner].getClusterIndex(order, idx).getIndex();}
    vector<vec> getPrimitiveClusterCoord(const primitive_position& pos) const { return _primitive_spin_clusters[pos.getCorner()].getClusterCoord( pos.getOrder(), pos.getIndex() );}
    int         getPrimitiveClusterNumber(int corner, int order) const { return _primitive_spin_clusters[corner].getClusterNum(order);}

    int         getGlobalPosition(int unit_cell_index, const primitive_position& pos);
    int         getGlobalClusterNumber(int order) const {return _lattice.getUnitCellNumber() * _primitive_cluster_size_fix_order(order);}
    uvec        getClusterIndex(int unit_cell, int corner_idx, int order, int idx);


private:
    void generate_primitive_clusters();

    cSpinCollection        _bath_spins;
    vector<cSPIN>          _spin_list;
    cDepthFirstPathTracing _primitive_dfpt;
    vector<cSpinCluster>   _primitive_spin_clusters;
    umat                   _primitive_cluster_size;
    urowvec                _primitive_cluster_size_fix_order;
    Lattice                _lattice;
    vector<int>            _center;
    int                    _unit_cell_num;
    int                    _atom_num_in_cell;

    primitive_position compute_sub_pos(const primitive_position& pos, size_t remove_k);
    vector<vec> relative_coordinates(const primitive_position& pos);
    pair< vector<int>, vector<vec> > remove_k(const primitive_position& pos, size_t remove_k);
    double patten_diff(const vector<vec>& v1, const vector<vec>& v2);
};
//}}}
////////////////////////////////////////////////////////////////////////////////
