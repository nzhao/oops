#include <armadillo>
#include "include/spin/SpinCollection.h"
#include "include/misc/misc.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinCollection
cSpinCollection::cSpinCollection()
{ //LOG(INFO) << "Default constructor of cSpinCollection.";
}
cSpinCollection::cSpinCollection(cSpinSource * source)
{
    _source = source;
}

cSpinCollection::~cSpinCollection()
{
    if (!_source) delete _source;
}

vector<cSPIN> cSpinCollection::getSpinList(const cClusterIndex& clst) const
{
    vector<cSPIN> sl;
    uvec id = clst.getIndex();
    for(int i=0; i<id.size(); ++i)
        sl.push_back(_spin_list[ id[i] ]);
    return sl;
}
void cSpinCollection::make()
{
/// call the 'generate' method of the cSpinSource to generate _spin_list.
    _spin_list= _source->generate();

    size_t nspin=_spin_list.size();
    mat d(nspin, nspin); d.zeros();
    for (int i=0; i<_spin_list.size(); ++i)
    {
        for (int j=i; j<_spin_list.size(); ++j)
        {
            //d(i,j)=distance(_spin_list[i].get_coordinate(),
            //                _spin_list[j].get_coordinate() );
            d(i,j)=spin_distance(_spin_list[i], _spin_list[j]);
        }
    }
    d =d+d.t();
    dist_mat=d;
}

sp_mat cSpinCollection::getConnectionMatrix(double threshold) const
{ 
    umat u = (dist_mat<=threshold);
    sp_mat res( conv_to< mat >::from(u) );
    for(int i=0; i<res.n_rows; ++i) res(i, i)=0;
    return res;
}

mat cSpinCollection::getCoordinateMat() const
{
    mat res=zeros<mat> (_spin_list.size(), 3);
    for(int i=0; i<_spin_list.size(); ++i)
        res.row(i) = trans( _spin_list[i].get_coordinate() );
    return res;
}


//}}}
////////////////////////////////////////////////////////////////////////////////
