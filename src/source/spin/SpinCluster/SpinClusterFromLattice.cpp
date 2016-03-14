#include "include/spin/SpinClusterFromLattice.h"

double PATTEN_EPS = 1e-10;

ostream& operator << (ostream& outs, const primitive_position& pos)
{
    outs << "corner = " << pos._corner << ", order = " << pos._order << ", index = " << pos._index;
    return outs;
}

cUniformBathOnLattice::cUniformBathOnLattice(const sp_mat& connection_matrix, size_t maxOrder, const cSpinCollection& bath_spins, const Lattice& lattice)
{
    _max_order = maxOrder;
    _connection_matrix=connection_matrix;
    _nspin     = connection_matrix.n_cols;
    _lattice   = lattice;
    _bath_spins = bath_spins;
    _spin_list = bath_spins.getSpinList();
    _has_cluster_index_list = false;
}

void cUniformBathOnLattice::generate()
{
    generate_primitive_clusters();

    //for(int corner=0; corner<_primitive_spin_clusters.size(); ++corner)
        //for(int order = 0; order<_max_order; ++order)
            //for(int idx = 0; idx < getPrimitiveClusterNumber(corner, order); ++idx)
            //{
                ////int corner = 0; int order =3; int idx = 2;
                //primitive_position pos(corner, order, idx);
                //cout << "pos= " << getPrimitiveClusterIndex(pos).t() ;
                //vector<primitive_position> pos_list = find_sub_position(pos);
                //for(int i=0; i< pos_list.size(); ++i)
                //{
                    //if(pos_list[i].getIndex() >=0)
                        //cout << "\t" << "sub_pos ="  << getPrimitiveClusterIndex(pos_list[i]).t();
                    //else
                    //{
                        //cout << "not found" <<endl;
                        //assert(0);
                    //}
                //}
            //}
}

void cUniformBathOnLattice::generate_primitive_clusters()
{
    vector<int> center = _lattice.getCenterSingleIndex();
    for(int i=0; i<center.size(); ++i)
    {
        mat init_mat = zeros<mat>(1, _nspin);
        init_mat(0, center[i]) = 1;
    
        _connection_matrix(span(0, center[i]-1), span::all) = zeros<mat>(center[i], _nspin);
        _connection_matrix(span::all, span(0, center[i]-1)) = zeros<mat>(_nspin, center[i]);

        cDepthFirstPathTracing dfpt_i(_connection_matrix, _max_order, init_mat);
        cSpinCluster spin_clusters_i(_bath_spins, &dfpt_i);
        spin_clusters_i.make();
        spin_clusters_i.diable_sub_cluster_position();
        cout << spin_clusters_i << endl;
        _primitive_spin_clusters.push_back(spin_clusters_i);
    }
    cout << "here" << endl;
    _primitive_cluster_size = zeros<umat>(center.size(), _max_order);
    for(int i=0; i<center.size(); ++i)
        for(int j=0; j<_max_order; ++j)
            _primitive_cluster_size(i, j) = getPrimitiveClusterNumber(i, j);
    cout << _primitive_cluster_size;
    _primitive_cluster_size_fix_order = sum(_primitive_cluster_size, 0);
    cout << _primitive_cluster_size_fix_order;

}

vector<primitive_position> cUniformBathOnLattice::find_sub_position(const primitive_position& pos)
{
    vector<primitive_position> pos_list;
    if(pos.getOrder() > 0)
    {
        for(int i=0; i<pos.getOrder()+1; ++i)
            pos_list.push_back( compute_sub_pos(pos, i) );
    }
    return pos_list;
}

primitive_position cUniformBathOnLattice::compute_sub_pos(const primitive_position& pos, size_t k)
{
    int corner = pos.getCorner(); int order = pos.getOrder(); int idx = pos.getIndex();
    pair< vector<int>, vector<vec> > p = remove_k(pos, k);
    vector<int> sub_idx = p.first;   vector<vec> sub_coord = p.second;

    int corner_idx = _lattice.getIndex( sub_idx[0]).back();
    int num_to_be_matched = getPrimitiveClusterNumber(corner_idx, order-1);

    int sub_pos = 1;
    for(int i=0; i<num_to_be_matched; ++i)
    {
        primitive_position pos_i(corner_idx, order-1, i);
        double diff = patten_diff( relative_coordinates(pos_i), sub_coord );
        if(diff < PATTEN_EPS)
        {
            sub_pos = i;
            break;
        }
    }

    primitive_position pos1( corner_idx,  order-1,  sub_pos);
    return pos1;
}


uvec cUniformBathOnLattice::getPrimitiveClusterIndex(const primitive_position& pos) const
{
    int corner = pos.getCorner(); int order = pos.getOrder(); int index = pos.getIndex();
    return getPrimitiveClusterIndex(corner, order, index);
}

vector<vec> cUniformBathOnLattice::relative_coordinates(const primitive_position& pos)
{
    vector<vec> res;
    vector<vec> clst_coord = getPrimitiveClusterCoord(pos);
    for(int i=0; i<clst_coord.size(); ++i)
        res.push_back( clst_coord[i] - clst_coord[0] );
    return res;
}

pair< vector<int>, vector<vec> >  cUniformBathOnLattice::remove_k(const primitive_position& pos, size_t k)
{
    uvec index = getPrimitiveClusterIndex(pos);
    vector<int> index_remove_k;
    for(int i=0; i<index.n_elem; ++i)
        if(i!=k)
            index_remove_k.push_back( index(i) );

    vec coord0 = _spin_list[index_remove_k[0]].get_coordinate();

    vector<vec> coord_remove_k;
    for(int i=0; i<index_remove_k.size(); ++i)
    {
        vec coord_i =  _spin_list[index_remove_k[i]].get_coordinate();
        coord_remove_k.push_back( coord_i - coord0 );
    }
    return make_pair( index_remove_k, coord_remove_k);
}

double cUniformBathOnLattice::patten_diff(const vector<vec>& v1, const vector<vec>& v2)
{
    double res = 0.0;
    assert( v1.size() == v2.size() );
    for(int i=0; i<v1.size(); ++i)
        res += norm( v1[i] - v2[i] );
    return res;
}

int cUniformBathOnLattice::global_position(int unit_cell_index, const primitive_position& pos)
{
    int order = pos.getOrder(); int corner = pos.getCorner(); int idx = pos.getIndex();
    int res = unit_cell_index * _primitive_cluster_size_fix_order(order)
        + corner * _primitive_cluster_size(corner, order) + idx;
    return res;
}
