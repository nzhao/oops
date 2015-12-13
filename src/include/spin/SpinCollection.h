#ifndef SPINCOLLECTION_H
#define SPINCOLLECTION_H

#include <vector>
#include <armadillo>
#include "include/spin/Spin.h"
#include "include/spin/SpinSource.h"

class cSpinCollection
{
public:
    cSpinCollection(cSpinSource * source);
    ~cSpinCollection();

    void make();
    vector<cSPIN>& getSpinList(){return spin_list;};
    arma::mat& getDistanceMatrix(){return dist_mat;};
    arma::umat getConnectionMatrix(double threshold);
private:
    cSpinSource* _source;
    vector<cSPIN> spin_list;

    arma::mat dist_mat;
};
#endif
