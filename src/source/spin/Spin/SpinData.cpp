#include <string>
#include "include/spin/SpinData.h"

SpinProperty C13      = {2,  6.728284e7,     0.0,      0.0};
SpinProperty N14      = {3,  1.9337792e7,    2.96,    0.0};
SpinProperty N15      = {2, -2.71261804e7,   0.0,      0.0};
SpinProperty ELECTRON = {2, -1.760859708e11, 0.0,      0.0};
SpinProperty NVE      = {3, -1.760859708e11, 2.87e3,   0.0};
SpinProperty Y        = {2, -1.3155e7,       0.0,      0.0};
SpinProperty Si       = {2, -5.3188e7,       0.0,      0.0};
////////////////////////////////////////////////////////////////////////////////
//{{{ SPINDATA
cSPINDATA::cSPINDATA()
{
    data["13C"]=C13;
    data["14C"]=N14;
    data["15C"]=N15;
    data["E"]=ELECTRON;
    data["NVe"]=NVE;
    data["Y"]=Y;
    data["Si"]=Si;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
