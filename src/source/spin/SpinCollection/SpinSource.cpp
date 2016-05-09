#include <vector>
#include <assert.h>
#include <string>
#include <fstream>
#include <iostream>
#include "include/spin/Spin.h"
#include "include/spin/SpinSource.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSource
cSpinSource::cSpinSource()
{ //LOG(INFO) << "Default constructor: cSpinSource.";
}

cSpinSource::~cSpinSource()
{ //LOG(INFO) << "Default destructor: cSpinSource.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSourceFromFile
cSpinSourceFromFile::cSpinSourceFromFile()
{ //LOG(INFO) << "Default constructor: cSpinSourceFromFile.";
}
cSpinSourceFromFile::cSpinSourceFromFile(string filename)
{ _filename=filename; }

cSpinSourceFromFile::~cSpinSourceFromFile()
{ //LOG(INFO) << "Default destructor: cSpinSourceFromFile.";
}

vector<cSPIN>& cSpinSourceFromFile::generate()
{
    int nbath;
    string atom_type;
    double x, y, z;

    ifstream coord(_filename.c_str());
    if(coord.fail())
     {
        cout<< "Input spin source opening failed."<<endl;
        if(!coord) assert(0);
     }
    

    coord >> nbath;
    for (int i = 0; i < nbath; ++i)
    {
        coord >> atom_type >> x >> y >> z; 
        vector<double> xyz; xyz.push_back(x); xyz.push_back(y); xyz.push_back(z);
        spin_list.push_back( cSPIN(xyz, atom_type) );
    }

    coord.close();
    return spin_list;
}
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSourceFromLattice
cSpinSourceFromLattice::cSpinSourceFromLattice(const Lattice& lattice,  const imat& range)
{
    _lattice = lattice;
    _lattice.setRange(range);
}

cSpinSourceFromLattice::cSpinSourceFromLattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const, int atom_num, const vector<vec>& pos, const vector<string>& isotope, const imat& range)
{
    _lattice = Lattice(dim, bases, lattice_const, atom_num, pos, isotope); 
    _lattice.setRange(range);
}

vector<cSPIN>& cSpinSourceFromLattice::generate()
{
    for(int i=0; i<_lattice.getTotalAtomNumber(); ++i)
    {
        cSPIN s( _lattice.getCoordinate(i), _lattice.getIsotope(i) );
        spin_list.push_back(s);
    }
    return spin_list;
}
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSourceUniformRandom
vector<cSPIN>& cSpinSourceUniformRandom::generate()
{
    arma_rng::set_seed(_seed);
    mat coord_mat(3, _max_num, fill::randu);

    vec v(_max_num);
    for(int i=0; i<_max_num; ++i)
        v[i]=norm(2.0*_range*coord_mat.col(i) - _range);
    uvec indices = sort_index(v);

    for(int i=0; i<_num; ++i)
    {
        vec coord = 2.0*_range*coord_mat.col( indices[i] ) - _range;
        cSPIN s(coord, _isotope );
        spin_list.push_back(s);
    }
    return spin_list;
}
//}}}
////////////////////////////////////////////////////////////////////////////////

