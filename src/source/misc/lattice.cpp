#include "include/misc/lattice.h"
////////////////////////////////////////////////////////////////////////////////
//{{{ Lattice
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

void Lattice::setRange(int max_idx)
{
    int absMax = max_idx > 0 ? max_idx : -max_idx;

    imat range(_dimension, 2);
    for(int i=0; i< _dimension; ++i)
    {
        range(i, 0) = -absMax;
        range(i, 1) = absMax + 1; 
    }
    setRange(range);
}

void Lattice::setRange(const imat& range)
{
    if(range.n_rows < _dimension)
        cout << "error range" <<endl;
    else
    {
        vector<int> temp_range_width;
        vector< pair<int, int> > temp_range;

        _unit_cell_num = 1;
        for(int i=0; i<_dimension; ++i)
        {
            int range_width_i = range(i,1) - range(i,0);
            _unit_cell_num *= range_width_i;
            temp_range_width.push_back( range_width_i );
            temp_range.push_back( make_pair(range(i, 0), range(i, 1) ) );
        }
        _total_atom_num = _atom_num_in_cell*_unit_cell_num;
        _range_width = temp_range_width;
        _range = temp_range;
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

int Lattice::getSingleIndex(const vector<int>& idx) const
{
    
    vector<int> base = _range_width;
    base.push_back( _atom_num_in_cell );
    vector<int> idx1;
    for(int i=0; i<_dimension; ++i)
        idx1.push_back( idx[i] - _range[i].first);
    idx1.push_back(idx[_dimension]);
    return base_number(idx1, base);
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
string Lattice::getIsotope(int i) const
{
    vector<int> idx = getIndex(i);
    return _isotope[ idx[_dimension] ];
}

vector< vector<int> > Lattice::getCenterIndex() const
{
    vector< vector<int> > res;
    for(int i=0; i<_atom_num_in_cell; ++i)
    {
        vector<int> res_i;
        for(int j=0; j<_dimension; ++j)
            res_i.push_back( _range[j].first + _range_width[j]/2 ); 
        res_i.push_back(i);
        res.push_back( res_i );
    }
    return res;
}

vector<int> Lattice::getCenterSingleIndex() const
{
    vector<int> res;
    vector< vector<int> > idx_list =  getCenterIndex();
    for(int i=0; i< idx_list.size(); ++i)
        res.push_back( getSingleIndex(idx_list[i]) );
    return res;
}

void Lattice::save_to_file(string filename)
{
    ofstream xyz(filename.c_str());
    if(!xyz) assert(0);

    xyz << _total_atom_num << endl;
    for(int i=0; i<_total_atom_num; ++i)
    {
        xyz << getIsotope(i) << "\t";
        vec coord = getCoordinate(i); 
        xyz << coord(0) << "\t" <<  coord(1) << "\t" << coord(2);
        xyz << endl;
    }
    xyz.close();

}

ostream&  operator << (ostream& outs, const Lattice& lattice)
{/*{{{*/
    outs << "dimension = " << lattice._dimension << endl;;
    for(int i=0; i<lattice._dimension; ++i)
        outs << "\t dim_" << i << " = [ " << lattice._range[i].first << ", " << lattice._range[i].second << " ), range_width = " << lattice._range_width[i]  << endl;

    outs<< "unit cell number = " << lattice._unit_cell_num << endl;
    outs<< "atom number per unit cell = " << lattice._atom_num_in_cell << endl;
    outs << "\t total atom number = " << lattice._total_atom_num << " = [ ";
    for(int i=0; i<lattice._dimension; ++i)
        outs << lattice._range_width[i] << " * ";
    outs << lattice._atom_num_in_cell << " ]" << endl;

    int num;
    if(lattice._total_atom_num > 20)
    {
        outs << "fisrt 20 atoms are listed below" << endl;
        num = 20;
    }
    else
    {
        outs << "all atoms are listed below" << endl;
        num = lattice._total_atom_num;
    }

    for(int i=0; i<num; ++i)
    {
        outs << "\t" << i << ": [ ";
        vector<int> idx = lattice.getIndex(i);
        for(int j=0; j<idx.size(); ++j)
        {
            outs<< idx[j] ; 
            if(j<idx.size() -1)
                outs << ", ";
        }
        outs << " ] = ";
        outs << lattice.getCoordinate(i).t();
        //outs << endl;
    }

    outs << "Lattice is centered at = " ;
    vector< vector<int> > idxC = lattice.getCenterIndex();
    for(int j=0; j<idxC.size(); ++j)
    {
        print_vector(idxC[j]);
        outs << " ";
    }
    outs <<  " = ";
    print_vector(lattice.getCenterSingleIndex() );
    return outs;
}/*}}}*/

//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ TwoDimFaceCenterLattice
TwoDimFaceCenterLattice::TwoDimFaceCenterLattice(double lattice_const, const vector<string>& isotope)
{
    _dimension = 2;
    vec base1, base2; base1 << 1.0 << 0.0 << 0.0; base2 << 0.0 << 1.0 << 0.0;
    _bases.push_back(base1); _bases.push_back(base2);
    _atom_num_in_cell = 2;
    vec coord1, coord2;
    coord1 << 0.0 << 0.0 << 0.0;
    coord2 << 0.5*lattice_const << 0.5*lattice_const<< 0.0;
    _pos_in_cell.push_back(coord1); _pos_in_cell.push_back(coord2);
    _lattice_constant.push_back(lattice_const); _lattice_constant.push_back(lattice_const);
    _isotope = isotope;
}

TwoDimFaceCenterLattice::TwoDimFaceCenterLattice(double lattice_const, const string& isotope)
{
    _dimension = 2;
    vec base1, base2; base1 << 1.0 << 0.0 << 0.0; base2 << 0.0 << 1.0 << 0.0;
    _bases.push_back(base1); _bases.push_back(base2);
    _atom_num_in_cell = 2;
    vec coord1, coord2;
    coord1 << 0.0 << 0.0 << 0.0;
    coord2 << 0.5*lattice_const << 0.5*lattice_const<< 0.0;
    _pos_in_cell.push_back(coord1); _pos_in_cell.push_back(coord2);
    _lattice_constant.push_back(lattice_const); _lattice_constant.push_back(lattice_const);
    for(int i=0; i<_atom_num_in_cell; ++i)
        _isotope.push_back(isotope);
}
//}}}
////////////////////////////////////////////////////////////////////////////////
