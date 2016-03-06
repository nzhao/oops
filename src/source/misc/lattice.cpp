#include "include/misc/lattice.h"

Lattice::Lattice(int dim, const vector<vec>& bases)
{
    _dimension = dim;
    _bases = bases;
    _atom_num_in_cell = 1;
    vec zero_vec = zeros<vec>(3); 
    _pos_in_cell.push_back( zero_vec );
}

Lattice::Lattice(int dim, const vector<vec>& bases, int atom_num, const vector<vec>& pos)
{
    _dimension = dim;
    _bases = bases;
    _atom_num_in_cell = atom_num;
    _pos_in_cell = pos;
}

void Lattice::setRange(const umat& range)
{
    if(range.n_rows < _dimension)
        cout << "error range" <<endl;
    else
    {
        for(int i=0; i<_dimension; ++i)
        {
            int range_width_i = range(i,1) - range(i,0);
            _range_width.push_back( range_width_i );
            _range.push_back( make_pair(range(i, 0), range(i, 1) ) );
        }
    }
}

vector<int> Lattice::getIndex(int num) const
{
    vector<int > res;
    res.push_back( num % _atom_num_in_cell );
    int q = num / _atom_num_in_cell;
    for(int i=0; i<_dimension; ++i)
    {
        res.push_back( _range[i].first + q % _range_base_num[i] );
        q = q / _range_base_num[i];
    }
    return res;
}
