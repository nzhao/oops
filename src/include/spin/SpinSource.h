#ifndef SPINSOURCE_H
#define SPINSOURCE_H
#include <vector>
#include <string>
#include "include/easylogging++.h"
#include "include/spin/Spin.h"

using namespace std;
/// \addtogroup SpinCollection
/// @{

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSource
/// This is an abstract class for generating spin_list
///
class cSpinSource
{
public:
    cSpinSource();
    virtual ~cSpinSource();
    virtual vector<cSPIN>& generate()=0;

    vector<cSPIN>& get_spin_list() {return spin_list;};
protected:
    vector<cSPIN> spin_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSourceFromFile
/// This class implements the spin_list generation from a given file.
/// 
class cSpinSourceFromFile:public cSpinSource
{
public:
    cSpinSourceFromFile();
    cSpinSourceFromFile(string filename);
    virtual ~cSpinSourceFromFile();

    vector<cSPIN>& generate();

private:
    void read_file();

    string _filename;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
#endif
