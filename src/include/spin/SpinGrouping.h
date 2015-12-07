#ifndef SPINGROUPING_H
#define SPINGROUPING_H
#include <vector>
#include <string>
#include "include/spin/Spin.h"

typedef vector<vector<int>> CLST_IDX;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSource
class cSpinGrouping
{
public:
    cSpinGrouping();
    virtual ~cSpinGrouping();
    virtual void generate()=0;

    vector<cSPIN>& get_cluster_index() {return cluster_index;};
protected:
    CLST_IDX cluster_index;
private:
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSourceFromFile
class cDepthFirstPathTracing:public cSpinGrouping
{
public:
    cDepthFirstPathTracing();
    cDepthFirstPathTracing(connection_matrix);
    virtual ~cDepthFirstPathTracing();

    void generate();

private:

};
#endif
