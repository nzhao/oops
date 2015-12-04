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
class cSpinGathering
{
public:
    cSpinGathering();
    virtual ~cSpinGathering();
    virtual void generate()=0;

    vector<cSPIN>& get_cluster_index() {return cluster_index;};
protected:
    CLST_IDX cluster_index;
private:
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSourceFromFile
class cSpinSourceFromFile:public cSpinSource
{
public:
    cSpinSourceFromFile();
    cSpinSourceFromFile(string filename);
    virtual ~cSpinSourceFromFile();

    void generate();

private:
    void read_file();

    string _filename;
};
#endif
