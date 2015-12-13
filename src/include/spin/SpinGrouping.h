#ifndef SPINGROUPING_H
#define SPINGROUPING_H
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <armadillo>
#include "include/spin/Spin.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cClusterIndex
class cClusterIndex
{
public:
    cClusterIndex();
    cClusterIndex(const vector<int>& idx);
    ~cClusterIndex();

    friend bool operator == (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend bool operator < (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend ostream&  operator << (ostream& outs, const cClusterIndex& idx);
private:
    vector<int> _index;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinGrouping

typedef set<cClusterIndex> CLST_IDX_LIST;

class cSpinGrouping
{
public:
    cSpinGrouping(arma::umat connection_matrix);
    cSpinGrouping();
    virtual ~cSpinGrouping();
    virtual CLST_IDX_LIST generate()=0;

    CLST_IDX_LIST get_cluster_index() {return _cluster_index_list;};
protected:
    CLST_IDX_LIST _cluster_index_list;
    arma::umat _connection_matrix;
private:
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSourceFromFile
class cDepthFirstPathTracing:public cSpinGrouping
{
public:
    cDepthFirstPathTracing();
    cDepthFirstPathTracing(arma::umat connection_matrix);
    virtual ~cDepthFirstPathTracing();

    CLST_IDX_LIST generate();

private:
    size_t max_size;

};
#endif
