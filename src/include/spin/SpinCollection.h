#ifndef SPINCOLLECTION_H
#define SPINCOLLECTION_H

#include <vector>
#include <armadillo>
#include "include/spin/Spin.h"
#include "include/spin/SpinSource.h"
//#include "include/spin/SpinClusterAlgorithm.h"
#include "include/spin/SpinClusterIndex.h"

using namespace std;
using namespace arma;
/// \addtogroup Spin
/// @{

/// \defgroup SpinList SpinList
/// @{

/// \defgroup SpinCollection SpinCollection
/// @{
////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinCollection
/// This class generates a collection of spins (spin_list) from a given cSpinSource.
/// This class also computes distances between spins, and connection matrix with a given threshlod distance.
///
class cSpinCollection
{
public:
    cSpinCollection();
    cSpinCollection(cSpinSource * source);
    ~cSpinCollection();

    //@{
    void make();
    //@}

    //@{
    size_t           getSpinNum() const {return _spin_list.size();};
    vector<cSPIN> getSpinList() const {return _spin_list;};
    vector<cSPIN> getSpinList(const cClusterIndex& clst) const;
    mat getCoordinateMat() const;
    mat& getDistanceMatrix() {return dist_mat;};
    sp_mat getConnectionMatrix (double threshold) const;
    //@}
private:
    cSpinSource* _source;
    vector<cSPIN> _spin_list;

    mat dist_mat;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}

/// @}
/// @}
#endif
