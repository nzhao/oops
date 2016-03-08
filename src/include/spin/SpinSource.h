#ifndef SPINSOURCE_H
#define SPINSOURCE_H
#include <vector>
#include <string>
#include "include/easylogging++.h"
#include "include/spin/Spin.h"
#include "include/misc/lattice.h"

using namespace std;
/// \addtogroup SpinCollection
/// @{

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSource
/// This is an abstract class for generating spin_list
///
class cSpinSource
{
public:
    cSpinSource();
    virtual ~cSpinSource();
    virtual vector<cSPIN>& generate()=0;

    vector<cSPIN>& get_spin_list() {return spin_list;};
protected:
    vector<cSPIN> spin_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSourceFromFile
/// This class implements the spin_list generation from a given file.
/// 
class cSpinSourceFromFile:public cSpinSource
{
public:
    cSpinSourceFromFile();
    cSpinSourceFromFile(string filename);
    virtual ~cSpinSourceFromFile();

    vector<cSPIN>& generate();

private:
    void read_file();

    string _filename;
};
//}}}

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSourceFromLattice
class cSpinSourceFromLattice:public cSpinSource
{
public:
    cSpinSourceFromLattice(){};
    cSpinSourceFromLattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const, int atom_num, const vector<vec>& pos, const vector<string>& isotope, const umat& range);
    ~cSpinSourceFromLattice(){};

    vector<cSPIN>& generate();
protected:
private:
    Lattice _lattice;
};
///}}} 
////////////////////////////////////////////////////////////////////////////////
/// @}
#endif
