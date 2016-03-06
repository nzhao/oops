#include <vector>
#include <armadillo>

using namespace arma;
using namespace std;

class Lattice
{
public:
    Lattice() {};
    Lattice(int dim, const vector<vec>& bases);
    Lattice(int dim, const vector<vec>& bases, int atom_num, const vector<vec>& pos);
    ~Lattice() {};

    vector<int> getIndex(int i) const;
    void setRange(const umat& range);
protected:
private:
    int                 _dimension;
    int                 _atom_num_in_cell;
    vector<vec>         _pos_in_cell;
    vector<vec>         _bases;
    vector< pair<int, int> > _range;
    vector<int>         _range_width;
    vector<int>         _range_base_num;
};
