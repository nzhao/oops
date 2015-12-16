#ifndef SPINCOLLECTION_H
#define SPINCOLLECTION_H

#include <vector>
#include <armadillo>
#include "include/spin/Spin.h"
#include "include/spin/SpinSource.h"

using namespace std;
using namespace arma;

class cSpinCollection
{
public:
    cSpinCollection(cSpinSource * source);
    ~cSpinCollection();

    void make();
    vector<cSPIN>& getSpinList(){return spin_list;};
    mat& getDistanceMatrix(){return dist_mat;};
    sp_mat getConnectionMatrix(double threshold);
private:
    cSpinSource* _source;
    vector<cSPIN> spin_list;

    mat dist_mat;
};
#endif
