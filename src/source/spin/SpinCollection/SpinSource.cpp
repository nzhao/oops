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
{ LOG(INFO) << "Default constructor: cSpinSource.";
}

cSpinSource::~cSpinSource()
{ LOG(INFO) << "Default destructor: cSpinSource.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSourceFromFile
cSpinSourceFromFile::cSpinSourceFromFile()
{ LOG(INFO) << "Default constructor: cSpinSourceFromFile.";
}
cSpinSourceFromFile::cSpinSourceFromFile(string filename)
{
    _filename=filename;
    cout << "filename=" << _filename << endl;
}

cSpinSourceFromFile::~cSpinSourceFromFile()
{ LOG(INFO) << "Default destructor: cSpinSourceFromFile.";
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
//}}}
////////////////////////////////////////////////////////////////////////////////
