#ifndef SPINCOLLECTION_H
#define SPINCOLLECTION_H

#include <vector>
#include "include/spin/Spin.h"
#include "include/spin/SpinSource.h"

class cSpinCollection
{
public:
    cSpinCollection(cSpinSource * source);
    ~cSpinCollection();

    void make();
    vector<cSPIN>& getSpinList(){return _source->get_spin_list();};
private:
    cSpinSource* _source;
};
#endif
