#ifndef SPINSOURCE_H
#define SPINSOURCE_H
#include <vector>
#include <string>
#include "include/spin/Spin.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinSource
class cSpinSource
{
public:
    cSpinSource();
    virtual ~cSpinSource();
    virtual void generate()=0;

    vector<cSPIN>& get_spin_list() {return spin_list;};
protected:
    vector<cSPIN> spin_list;
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
