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
//{{{
cSpinSourceFromLattice::cSpinSourceFromLattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const, int atom_num, const vector<vec>& pos, const vector<string>& isotope, const umat& range)
{
    _lattice = Lattice(dim, bases, lattice_const, atom_num, pos, isotope); 
    _lattice.setRange(range);
    _lattice.generate_spins();
}

vector<cSPIN>& cSpinSourceFromLattice::generate()
{
    spin_list = _lattice.getSpinList();
    return spin_list;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
