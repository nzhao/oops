#include "include/misc/lattice.h"

Lattice::Lattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const)
{
    _dimension = dim;
    _bases = bases;
    _atom_num_in_cell = 1;
    vec zero_vec = zeros<vec>(3); 
    _pos_in_cell.push_back( zero_vec );
    _lattice_constant = lattice_const;
}

Lattice::Lattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const, int atom_num, const vector<vec>& pos, const vector<string>& isotope)
{
    _dimension = dim;
    _bases = bases;
    _atom_num_in_cell = atom_num;
    _pos_in_cell = pos;
    _lattice_constant = lattice_const;
    _isotope = isotope;
}

void Lattice::setRange(const umat& range)
{
    if(range.n_rows < _dimension)
        cout << "error range" <<endl;
    else
    {
        _total_atom_num = _atom_num_in_cell;
        for(int i=0; i<_dimension; ++i)
        {
            int range_width_i = range(i,1) - range(i,0);
            _total_atom_num *= range_width_i;
            _range_width.push_back( range_width_i );
            _range.push_back( make_pair(range(i, 0), range(i, 1) ) );
        }
    }
}

vector<int> Lattice::getIndex(int num) const
{
    vector<int> base = _range_width;
    base.push_back( _atom_num_in_cell );
    vector<int> res = base_transform(num, base);
    for(int i=0; i<_dimension; ++i)
        res[i] += _range[i].first;
    return res;
}

vec Lattice::getCoordinate(const vector<int>& idx) const
{
    vec coord = zeros<vec>(3);
    for(int i=0; i<_dimension; ++i)
        coord += idx[i] * _lattice_constant[i] * _bases[i];
    coord += _pos_in_cell[ idx[_dimension] ];
    return coord;
}
vec Lattice::getCoordinate(int i) const
{
    vector<int> idx = getIndex(i);
    vec coord = getCoordinate(idx);
    return coord; 
}

void Lattice::generate_spins()
{
    for(int i=0; i<_total_atom_num; ++i)
    {
        vector<int> idx = getIndex(i);
        int idx_in_cell = idx[_dimension];
        cSPIN s(getCoordinate(idx), _isotope[idx_in_cell]);
        _spin_list.push_back(s);
    }
}
