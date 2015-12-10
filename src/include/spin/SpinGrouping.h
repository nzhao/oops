#ifndef SPINGROUPING_H
#define SPINGROUPING_H
#include <vector>
#include <string>
#include "include/spin/Spin.h"

typedef vector<vector<int>> VECT2;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSource
class cSpinGrouping
{
public:
    cSpinGrouping(VECT2 connection_matrix);
    virtual ~cSpinGrouping();
    virtual void generate()=0;

    VECT2& get_cluster_index() {return cluster_index;};
protected:
    VECT2 _cluster_index;
    VECT2 _connection_matrix;
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
    size_t max_size;

};
#endif
