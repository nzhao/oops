#include "include/spin/SpinCollection.h"

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
    _source->generate();
}

