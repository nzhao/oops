#ifndef SPINCOLLECTION_H
#define SPINCOLLECTION_H

#include <vector>
#include "include/spin/Spin.h"
#include "include/spin/SpinSource.h"

class cSpinCollection
{
public:
    cSpinCollection(cSpinSource & source);
private:
    vector<cSPIN> spin_list;
    cSpinSource& spin_source;
};
#endif
