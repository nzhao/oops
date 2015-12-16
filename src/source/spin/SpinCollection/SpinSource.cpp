#include <vector>
#include <assert.h>
#include <string>
#include <fstream>
#include <iostream>
#include "include/spin/Spin.h"
#include "include/spin/SpinSource.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSource

cSpinSource::cSpinSource()
{
    vector<cSPIN> spin_list;
//    cout << "spin_list is initialized in cSpinSource" << endl;
}

cSpinSource::~cSpinSource()
{// cout << "destructor: cSpinSource is called" << endl;
}

vector<cSPIN>& cSpinSource::generate()
{ //cout << "cSpinSource::generate is called" << endl;
    return spin_list;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSourceFromFile

cSpinSourceFromFile::cSpinSourceFromFile()
{
    cout << "Need a file name!" << endl;
}
cSpinSourceFromFile::cSpinSourceFromFile(string filename)
{
    _filename=filename;
    cout << "filename=" << _filename << endl;
}

cSpinSourceFromFile::~cSpinSourceFromFile()
{// cout << "destructor: cSpinSourceFromFile is called" << endl;
}

vector<cSPIN>& cSpinSourceFromFile::generate()
{
    int nbath;
    string atom_type;
    double x, y, z;

    ifstream coord(_filename);
    if(!coord) assert(0);

    coord >> nbath;
    for (int i = 0; i < nbath; ++i)
    {
        coord >> atom_type >> x >> y >> z; 
        spin_list.push_back( cSPIN(vector<double> {x, y, z}, atom_type) );
    }

    coord.close();
    return spin_list;
}

