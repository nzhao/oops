#include <vector>
#include <armadillo>
#include "include/misc/misc.h"

using namespace arma;
using namespace std;

class Lattice
{
public:
    Lattice() {};
    Lattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const);
    Lattice(int dim, const vector<vec>& bases, const vector<double>& lattice_const, int atom_num, const vector<vec>& pos);
    ~Lattice() {};

    vector<int> getIndex(int i) const;
    vec         getCoordinate(int i) const;
    vec         getCoordinate(const vector<int>& idx) const;
    void        setRange(const umat& range);
protected:
private:
    int                 _dimension;
    int                 _atom_num_in_cell;
    vector<vec>         _pos_in_cell;
    vector<vec>         _bases;
    vector< pair<int, int> > _range;
    vector<int>         _range_width;
    vector<double>      _lattice_constant;
};
