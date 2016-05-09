#ifndef SPINSOURCE_H
#define SPINSOURCE_H
#include <vector>
#include <string>
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
    cSpinSourceFromLattice(const Lattice& lattice) {_lattice = lattice;};
    cSpinSourceFromLattice(const Lattice& lattice,  const imat& range);
    cSpinSourceFromLattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const, int atom_num, const vector<vec>& pos, const vector<string>& isotope, const imat& range);
    ~cSpinSourceFromLattice(){};

    vector<cSPIN>& generate();
protected:
private:
    Lattice _lattice;
};
///}}} 
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSourceUniformRandom
class cSpinSourceUniformRandom:public cSpinSource
{
public:
    cSpinSourceUniformRandom(const double range, const int max_num, const int num, const string& isotope, const int seed) {_range = range; _max_num = max_num; _num = num; _isotope = isotope; _seed = seed;}
    vector<cSPIN>& generate();
protected:
private:
    int _num;
    int _max_num;
    double _range;
    string _isotope;
    int _seed;
};
/// @}
#endif
