#ifndef LATTICE_H
#define LATTICE_H
#include <vector>
#include <armadillo>
#include "include/misc/misc.h"
#include "include/spin/Spin.h"

using namespace arma;
using namespace std;

class Lattice
{
public:
    Lattice() {};
    Lattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const);
    Lattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const, int atom_num, const vector<vec>& pos, const vector<string>& isotope);
    ~Lattice() {};

    vector<int> getIndex(int i) const;
    vec         getCoordinate(int i) const;
    vec         getCoordinate(const vector<int>& idx) const;
    vector<cSPIN> getSpinList() const {return _spin_list;};

    void        setRange(const umat& range);
    void        generate_spins();
    int         getSpinNum() const {return _spin_list.size();};
protected:
private:
    int                 _dimension;
    int                 _atom_num_in_cell;
    int                 _total_atom_num;
    vector<vec>         _pos_in_cell;
    vector<vec>         _bases;
    vector< pair<int, int> > _range;
    vector<int>         _range_width;
    vector<double>      _lattice_constant;
    vector<string>      _isotope;
    vector<cSPIN>       _spin_list;
};
#endif
