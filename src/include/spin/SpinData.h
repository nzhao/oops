#ifndef SPINDATA_H
#define SPINDATA_H

#include <string>
#include <map>

using namespace std;

/// \addtogroup Spin
/// @{

////////////////////////////////////////////////////////////////////////////////
//{{{ SpinProperty
/// Spin data.
///
struct SpinProperty
{
    int multiplicity; ///< 2S+1.
    double gamma; ///< in SI unit: rad/s/T.
    double omegaQ; ///< in MHz.
    double eta; ///< in MHz.
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSPINDATA
/// This class provides spin data.
/// The data are stored in a form defined by structure SpinProperty.
///
class cSPINDATA
{
public:
    cSPINDATA();
    SpinProperty getData(string name) {return data[name]; };
private:
    map<string, SpinProperty> data;

};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
#endif
