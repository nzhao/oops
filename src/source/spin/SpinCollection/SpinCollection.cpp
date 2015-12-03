#include "include/spin/SpinCollection.h"

cSpinCollection::cSpinCollection(cSpinSource & source)
{
    spin_source = source;
}

cSpinCollection::make()
{
    spin_source.generate();
}

