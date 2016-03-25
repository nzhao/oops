#include "include/spin/SpinClusterFromLattice.h"

double PATTEN_EPS = 1e-10;

cUniformBathOnLattice::cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins, const Lattice& lattice, int root_range_idx)
{
    _max_order         = maxOrder;
    _bath_spins        = bath_spins;
    _connection_matrix = connection_matrix;
    _nspin             = connection_matrix.n_cols;
    _spin_list         = bath_spins.getSpinList();

    _lattice          = lattice;
    _center           = lattice.getCenterSingleIndex();
    _atom_num_in_cell = lattice.getUnitCellAtomNumber();
    
    _root_range = zeros<imat> (_lattice.getDimension(), 2);
    for(int i=0; i< _lattice.getDimension(); ++i)
    {
        _root_range(i, 0) = -root_range_idx;
        _root_range(i, 1) = root_range_idx + 1; 
    }

    _root_lattice = _lattice;
    _root_lattice.setRange(_root_range); 
    _root_center = _root_lattice.getCenterSingleIndex();
    _unit_cell_num = _root_lattice.getUnitCellNumber();
    cout << _root_lattice << endl;
}

cUniformBathOnLattice::cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins, const Lattice& lattice, const imat& root_range)
{
    _max_order         = maxOrder;
    _bath_spins        = bath_spins;
    _connection_matrix = connection_matrix;
    _nspin             = connection_matrix.n_cols;
    _spin_list         = bath_spins.getSpinList();

    _lattice          = lattice;
    _center           = lattice.getCenterSingleIndex();
    _atom_num_in_cell = lattice.getUnitCellAtomNumber();
    _root_range       = root_range;

    _root_lattice = _lattice;
    _root_lattice.setRange(_root_range); 
    _root_center = _root_lattice.getCenterSingleIndex();
    _unit_cell_num = _root_lattice.getUnitCellNumber();
    cout << _root_lattice << endl;
}

void cUniformBathOnLattice::generate()
{
    generate_primitive_clusters();
    generate_sub_primitive_position();
    generate_cluster_index_list();

}

void cUniformBathOnLattice::generate_primitive_clusters()
{/*{{{*/
    for(int i=0; i<_atom_num_in_cell; ++i)
    {
        cout << "generating primitive clusters of atom " << i << "/" << _atom_num_in_cell << " ..." << endl;
        mat init_mat = zeros<mat>(1, _nspin);
        init_mat(0, _center[i]) = 1;
        
        _connection_matrix(span(0, _center[i]-1), span::all) = zeros<mat>(_center[i], _nspin);
        _connection_matrix(span::all, span(0, _center[i]-1)) = zeros<mat>(_nspin, _center[i]);
        
        cDepthFirstPathTracing dfpt_i(_connection_matrix, _max_order, init_mat);
        cSpinCluster spin_clusters_i(_bath_spins, &dfpt_i);
        spin_clusters_i.make();
        spin_clusters_i.diable_sub_cluster_position();
        //cout << spin_clusters_i << endl;
        _primitive_spin_clusters.push_back(spin_clusters_i);
    }

    for(int order=0; order<_max_order; ++order)
    {
        umat m_order = _primitive_spin_clusters[0].getClusterIndex(order); 
        for(int i=1; i<_atom_num_in_cell; ++i)
            m_order = join_vert(m_order, _primitive_spin_clusters[i].getClusterIndex(order) );
        _primitive_cluster_mat.push_back( m_order);
        //cout << m_order << endl;
    }

    _primitive_cluster_size = zeros<umat>(_atom_num_in_cell, _max_order);
    for(int i=0; i<_atom_num_in_cell; ++i)
        for(int j=0; j<_max_order; ++j)
            _primitive_cluster_size(i, j) = getPrimitiveClusterNumber(i, j);
    _primitive_cumsum_size = zeros<umat>(_atom_num_in_cell+1, _max_order);
    _primitive_cumsum_size(span(1,  _atom_num_in_cell), span::all) = cumsum(_primitive_cluster_size);

    for(int order=0; order<_max_order; ++order)
        _total_cluster_number.push_back( _unit_cell_num*_primitive_cumsum_size(_atom_num_in_cell, order) );
    cout << "Summary: cluster numbers" << endl << _primitive_cumsum_size << endl;

}/*}}}*/

void cUniformBathOnLattice::generate_sub_primitive_position()
{/*{{{*/
    for(int order=1; order<_max_order; ++order)
    {
        cout << "generating sub_primitive_position of order = " << order << "/" << _max_order << " ... " << endl;
        SubPosLst_FixOrder sp_fix_order;
        for(int i=0; i<_primitive_cluster_mat[order].n_rows; ++i)
        {
            urowvec v = _primitive_cluster_mat[order].row(i);
            //cout << "###############" << endl << "v = " << v;
            
            SubPosLst pos_list;
            for(int rm_j=0; rm_j<v.n_elem; ++rm_j)
            {
                urowvec v_rm_j = v; v_rm_j.shed_cols(rm_j, rm_j);
                //cout << "sub_v(" << rm_j << ") = " << v_rm_j;
                int pos = locate_primitive_sub_clusters(v_rm_j); 
                //cout << "pos = " << pos << endl;
                if(pos>=0)
                    pos_list.push_back( (size_t) pos );
            }
            //print_vector(pos_list);
            //cout << endl << endl;
            sp_fix_order.push_back( pos_list );
        }
        _sub_pos.push_back( sp_fix_order );
    }
}/*}}}*/

int cUniformBathOnLattice::locate_primitive_sub_clusters(const urowvec& v)
{/*{{{*/
    vector<int> global_lattice_idx = _lattice.getIndex( v(0) );

    int order = v.n_elem-1;
    int idx_in_unit_cell = global_lattice_idx.back();
    int root_single_idx = _root_lattice.getSingleIndex(global_lattice_idx);

    int cell_diff = (root_single_idx - _root_center[idx_in_unit_cell]) / _atom_num_in_cell;
    urowvec v_shift = v - v(0) + _center[idx_in_unit_cell]; 
    int i0 = _primitive_cumsum_size(idx_in_unit_cell, order);
    int i1 = _primitive_cumsum_size(idx_in_unit_cell+1, order);

    int res = cell_diff * _primitive_cumsum_size( _atom_num_in_cell, order);
    for(int i=i0; i<i1; ++i)
        if( all( v_shift==_primitive_cluster_mat[order].row(i) ) )
        {
            return res + i;
        }
    return -1;
}/*}}}*/

void cUniformBathOnLattice::generate_cluster_index_list()
{
    _cluster_index_list.clear();
    for(int order=0; order<_max_order; ++order)
    {
        cout << "generating cluster index list of order " << order << "/" << _max_order << " ... " << endl;
        FIX_ORDER_INDEX_SET sub_pos_set;

        int pr_clst_num = _primitive_cumsum_size(_atom_num_in_cell, order);
        int clst_num = pr_clst_num*_unit_cell_num;
        umat clst_mat = zeros<umat>(clst_num, order+1);
        for(int cell_idx=0; cell_idx<_unit_cell_num; ++cell_idx)
        {
            vector<int> vIdx = _root_lattice.getIndex(cell_idx*_atom_num_in_cell);
            int cell_shift = _lattice.getSingleIndex(vIdx) - _center[0];
            //for(int i=0; i<_primitive_cluster_mat[order].n_rows; ++i)
            for(int i=0; i<pr_clst_num; ++i)
            {
               urowvec idx = _primitive_cluster_mat[order].row(i) + cell_shift;
               cClusterIndex cIdx( idx.t() );
               clst_mat.row(cell_idx*pr_clst_num+i) = idx;
               if(order > 0)
               {
                   vector<size_t> shift_sub_pos;
                   for(int q=0; q<_sub_pos[order-1][i].size(); ++q)
                   {
                       //cout << _sub_pos[order-1][i][q] << ", " << cell_idx << endl;
                       int pos = _sub_pos[order-1][i][q] + cell_idx*_primitive_cumsum_size(_atom_num_in_cell, order-1);
                       //cout << "ppos = " << pos << endl;
                       //cout << "tot = " << _total_cluster_number[order-1] << endl;
                       if(pos < _total_cluster_number[order-1])
                           shift_sub_pos.push_back( pos ); 
                   }
                   //cIdx.setSubClstPos( _sub_pos[order-1][i] );
                   cIdx.setSubClstPos( shift_sub_pos );
               }
               sub_pos_set.insert( cIdx );
               //cout << cIdx << endl;
            }
        }
        _cluster_index_list.push_back( sub_pos_set );
        _cluster_index_mat.push_back(clst_mat);
        //cout << clst_mat.head_rows(100);
    }
}

