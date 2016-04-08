#ifndef LATTICE_H
#define LATTICE_H
#include <vector>
#include <armadillo>
#include "include/misc/misc.h"
#include "include/spin/Spin.h"

using namespace arma;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
//{{{ Lattice
class Lattice
{
public:
    Lattice() {};
    Lattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const);
    Lattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const, int atom_num, const vector<vec>& pos, const vector<string>& isotope);
    ~Lattice() {};

    int         getDimension() const {return _dimension;};
    vector<int> getIndex(int i) const;
    vec         getCoordinate(int i) const;
    vec         getCoordinate(const vector<int>& idx) const;
    vector<vec> getBases() const {return _bases;};
    vector<int> getRangeWidth() const {return _range_width;};
    string      getIsotope(int i) const;
    int         getTotalAtomNumber() const {return _total_atom_num;};
    int         getUnitCellNumber() const {return _unit_cell_num;};
    int         getUnitCellAtomNumber() const {return _atom_num_in_cell;};
    int         getSingleIndex(const vector<int>& idx) const;
    vector< vector<int> > getCenterIndex() const;
    vector<int> getCenterSingleIndex() const;

    void        setRange(int i);
    void        setRange(const imat& range);
    void        save_to_file(string filename);

    friend ostream&  operator << (ostream& outs, const Lattice& lattice);

protected:
    int                 _dimension;
    int                 _atom_num_in_cell;
    int                 _unit_cell_num;
    int                 _total_atom_num;
    vector<vec>         _pos_in_cell;
    vector<vec>         _bases;
    vector< pair<int, int> > _range;
    vector<int>         _range_width;
    vector<double>      _lattice_constant;
    vector<string>      _isotope;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ TwoDimFaceCenterLattice
class TwoDimFaceCenterLattice:public Lattice
{
public:
    TwoDimFaceCenterLattice(){};
    TwoDimFaceCenterLattice(double lattice_const, const vector<string>& isotops);
    TwoDimFaceCenterLattice(double lattice_const, const string& isotope);
    ~TwoDimFaceCenterLattice(){};
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
