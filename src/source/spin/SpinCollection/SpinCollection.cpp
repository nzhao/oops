#include <armadillo>
#include "include/spin/SpinCollection.h"
#include "include/misc/misc.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinCollection
cSpinCollection::cSpinCollection(cSpinSource * source)
{
    _source = source;
}

cSpinCollection::~cSpinCollection()
{
    if (!_source) delete _source;
}

void cSpinCollection::make()
{
/// call the 'generate' method of the cSpinSource to generate spin_list.
    spin_list= _source->generate();

    int nspin=spin_list.size();
    mat d(nspin, nspin); d.zeros();
    for (int i=0; i<spin_list.size(); ++i)
    {
        for (int j=i; j<spin_list.size(); ++j)
        {
            //d(i,j)=distance(spin_list[i].get_coordinate(),
            //                spin_list[j].get_coordinate() );
            d(i,j)=spin_distance(spin_list[i], spin_list[j]);
        }
    }
    d =d+d.t();
    dist_mat=d;
}

sp_mat cSpinCollection::getConnectionMatrix(double threshold)
{ 
    umat u = (dist_mat<=threshold);
    sp_mat res( conv_to< mat >::from(u) );
    for(int i=0; i<res.n_rows; ++i) res(i, i)=0;
    return res;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
