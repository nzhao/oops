#include "armadillo"
#include "include/spin/SpinCollection.h"
#include "include/misc/misc.h"

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
    spin_list= _source->generate();

    int nspin=spin_list.size();
    arma::mat d(nspin, nspin);
    for (int i=0; i<spin_list.size(); ++i)
    {
        for (int j=i; j<spin_list.size(); ++j)
        {
            d(i,j)=distance(spin_list[i],spin_list[j]);
        }
    }
    d =d+d.t();
    dist_mat=d;
}

arma::umat cSpinCollection::getConnectionMatrix(double threshold)
{
    return (dist_mat<=threshold);
}

