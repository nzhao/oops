#include "include/spin/SpinClusterFromLattice.h"

cUniformBathOnLattice::cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins)
{
    _max_order = maxOrder;
    _connection_matrix=connection_matrix;
    _nspin     = connection_matrix.n_cols;
    _bath_spins = bath_spins;
    _bath_center_index = _nspin/2;
}

void cUniformBathOnLattice::generate()
{
    generate_primitive_clusters();
}

void cUniformBathOnLattice::generate_primitive_clusters()
{

    mat init_mat = zeros<mat>(1, _nspin);
    init_mat(_bath_center_index) = 1;

    _connection_matrix(span(0, _bath_center_index-1), span::all) = zeros<mat>(_bath_center_index, _nspin);
    _connection_matrix(span::all, span(0, _bath_center_index-1)) = zeros<mat>(_nspin, _bath_center_index);

    _primitive_dfpt = cDepthFirstPathTracing(_connection_matrix, _max_order, init_mat);
    _primitive_spin_clusters = cSpinCluster(_bath_spins, &_primitive_dfpt);
    _primitive_spin_clusters.make();
    cout << _primitive_spin_clusters << endl;
}
