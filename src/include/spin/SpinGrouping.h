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
    cClusterIndex(const arma::uvec& idx);
    ~cClusterIndex();

    arma::sp_mat get_array(size_t nspin);

    friend bool operator == (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend bool operator < (const cClusterIndex& idx1, const cClusterIndex& idx2);
    friend ostream&  operator << (ostream& outs, const cClusterIndex& idx);
private:
    arma::uvec _index;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinGrouping

typedef set<cClusterIndex> CLST_IDX_LIST;
typedef vector<CLST_IDX_LIST> CLST;

class cSpinGrouping
{
public:
    cSpinGrouping(arma::umat connection_matrix);
    cSpinGrouping();
    virtual ~cSpinGrouping();
    virtual CLST generate()=0;

    CLST get_cluster_index() {return _cluster_index_list;};
    arma::sp_mat get_cluster_mat(int order);
protected:
    size_t nspin;
    arma::umat _connection_matrix;
    CLST _cluster_index_list;
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

    CLST generate();

private:
    size_t max_size;

};
#endif
