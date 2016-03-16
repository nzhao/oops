#include "include/spin/SpinClusterFromLattice.h"

double PATTEN_EPS = 1e-10;

cUniformBathOnLattice::cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins, const Lattice& lattice)
{
    _max_order         = maxOrder;
    _bath_spins        = bath_spins;
    _connection_matrix = connection_matrix;
    _nspin             = connection_matrix.n_cols;
    _spin_list         = bath_spins.getSpinList();

    _lattice          = lattice;
    _center           = lattice.getCenterSingleIndex();
    _atom_num_in_cell = lattice.getUnitCellAtomNumber();
    _unit_cell_num    = lattice.getUnitCellNumber();
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
{
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
                if(pos>0)
                    pos_list.push_back( (size_t) pos );
            }
            //print_vector(pos_list);
            //cout << endl << endl;
            sp_fix_order.push_back( pos_list );
        }
        _sub_pos.push_back( sp_fix_order );
    }
}

int cUniformBathOnLattice::locate_primitive_sub_clusters(const urowvec& v)
{
    int order = v.n_elem-1;
    int idx_in_unit_cell = _lattice.getIndex( v(0) ).back();
    int unit_cell_pos    = ( v(0) - idx_in_unit_cell ) / _atom_num_in_cell;
    int res              = unit_cell_pos * _primitive_cumsum_size( _atom_num_in_cell, order);

    urowvec v_shift = v - v(0) + _center[idx_in_unit_cell];
    int i0 = _primitive_cumsum_size(idx_in_unit_cell, order);
    int i1 = _primitive_cumsum_size(idx_in_unit_cell+1, order);

    //cout << "unit_cell_pos = " << unit_cell_pos << ", v_shift = " << v_shift;
    for(int i=i0; i<i1; ++i)
        if( all( v_shift==_primitive_cluster_mat[order].row(i) ) )
        {
            //cout << "( " << res << ", " << i << ")" << endl << endl;;
            return res + i;
        }
    //cout << "not found" << endl << endl;
    return -1;
}
void cUniformBathOnLattice::generate_cluster_index_list()
{
    _cluster_index_list.clear();
    for(int order=0; order<_max_order; ++order)
    {
        cout << "generating cluster index list of order " << order << "/" << _max_order << " ... " << endl;
        FIX_ORDER_INDEX_SET sub_pos_set;
        for(int cell_idx=0; cell_idx<_unit_cell_num; ++cell_idx)
        {
            int cell_shift = cell_idx*_atom_num_in_cell - _center[0];
            for(int i=0; i<_primitive_cluster_mat[order].n_rows; ++i)
            {
               urowvec idx = _primitive_cluster_mat[order].row(i) + cell_shift;
               cClusterIndex cIdx( idx.t() );
               if(order > 0)
               {
                   vector<size_t> shift_sub_pos;
                   for(int q=0; q<_sub_pos[order-1][i].size(); ++q)
                   {
                       int pos = _sub_pos[order-1][i][q]+ cell_shift;
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
    }

}

///{{{
//void cUniformBathOnLattice::generate_cluster_index_list()
//{
    //for(int order=0; order<_max_order; ++order)
    //{
        //cout << "order = " << order << " has " << getGlobalClusterNumber(order) << " clusters." << endl;
        //for(int i=0; i<_unit_cell_num; ++i)
        //{
            //for(int j=0; j<_atom_num_in_cell; ++j)
            //{
                //for(int k=0; k<_primitive_cluster_size(j, order); ++k)
                //{
                    //cClusterIndex cIdx( getClusterIndex(i, j, order, k) );
                    //_cluster_index_list[order].insert(cIdx);
                //}
            //}
        //}
    //}
//}
//uvec cUniformBathOnLattice::getClusterIndex(int unit_cell, int corner_idx, int order, int idx)
//{
    //uvec pr_v = getPrimitiveClusterIndex(corner_idx, order, idx);
    //uvec v = pr_v - _center[0] + unit_cell*_atom_num_in_cell;
    //return v;
//}

//vector<primitive_position> cUniformBathOnLattice::find_sub_position(const primitive_position& pos)
//{
    //vector<primitive_position> pos_list;
    //if(pos.getOrder() > 0)
        //for(int i=0; i<pos.getOrder()+1; ++i)
            //pos_list.push_back( compute_sub_pos(pos, i) );
    //return pos_list;
//}
//pair<int, int>  cUniformBathOnLattice::match_pos(const urowvec& idx)
//{
    //int order = idx.n_elem;
    //int vertex, pos;

    //int corner_idx = _lattice.getIndex( idx(0) ).back();
    //umat m = _primitive_spin_clusters[corner_idx].getClusterIndex(order);
    //for(int i=0; i< m.n_rows; ++i)
    //{
        //urowvec row_i = m.row(i) - m(i, 0);
        //urowvec idx0  = idx - idx(0);
        
    //}
    //return make_pair( vertex, pos);
//}

//primitive_position cUniformBathOnLattice::compute_sub_pos(const primitive_position& pos, size_t k)
//{
    //int corner = pos.getCorner(); int order = pos.getOrder(); int idx = pos.getIndex();
    //pair< vector<int>, vector<vec> > p = remove_k(pos, k);
    //vector<int> sub_idx = p.first;   vector<vec> sub_coord = p.second;

    //int corner_idx = _lattice.getIndex( sub_idx[0]).back();
    //int num_to_be_matched = getPrimitiveClusterNumber(corner_idx, order-1);

    //int sub_pos = 1;
    //for(int i=0; i<num_to_be_matched; ++i)
    //{
        //primitive_position pos_i(corner_idx, order-1, i);
        //double diff = patten_diff( relative_coordinates(pos_i), sub_coord );
        //if(diff < PATTEN_EPS)
        //{
            //sub_pos = i;
            //break;
        //}
    //}

    //primitive_position pos1( corner_idx,  order-1,  sub_pos);
    //return pos1;
//}

//uvec cUniformBathOnLattice::getPrimitiveClusterIndex(const primitive_position& pos) const
//{
    //int corner = pos.getCorner(); int order = pos.getOrder(); int index = pos.getIndex();
    //return getPrimitiveClusterIndex(corner, order, index);
//}
//vector<vec> cUniformBathOnLattice::relative_coordinates(const primitive_position& pos)
//{
    //vector<vec> res;
    //vector<vec> clst_coord = getPrimitiveClusterCoord(pos);
    //for(int i=0; i<clst_coord.size(); ++i)
        //res.push_back( clst_coord[i] - clst_coord[0] );
    //return res;
//}

//pair< vector<int>, vector<vec> >  cUniformBathOnLattice::remove_k(const primitive_position& pos, size_t k)
//{
    //uvec index = getPrimitiveClusterIndex(pos);
    //vector<int> index_remove_k;
    //for(int i=0; i<index.n_elem; ++i)
        //if(i!=k)
            //index_remove_k.push_back( index(i) );

    //vec coord0 = _spin_list[index_remove_k[0]].get_coordinate();

    //vector<vec> coord_remove_k;
    //for(int i=0; i<index_remove_k.size(); ++i)
    //{
        //vec coord_i =  _spin_list[index_remove_k[i]].get_coordinate();
        //coord_remove_k.push_back( coord_i - coord0 );
    //}
    //return make_pair( index_remove_k, coord_remove_k);
//}

//double cUniformBathOnLattice::patten_diff(const vector<vec>& v1, const vector<vec>& v2)
//{
    //double res = 0.0;
    //assert( v1.size() == v2.size() );
    //for(int i=0; i<v1.size(); ++i)
        //res += norm( v1[i] - v2[i] );
    //return res;
//}
//int cUniformBathOnLattice::getGlobalPosition(int unit_cell_index, const primitive_position& pos)
//{
    //int order = pos.getOrder(); int corner = pos.getCorner(); int idx = pos.getIndex();
    //int res = unit_cell_index * _primitive_cumsum_size(_atom_num_in_cell, order)
        //+ corner * _primitive_cluster_size(corner, order) + idx;
    //return res;
//}
//}}}

