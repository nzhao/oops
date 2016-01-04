#ifndef SPINDATA_H
#define SPINDATA_H

#include <string>
#include <map>

using namespace std;

struct SpinProperty
{
    int multiplicity;
    double gamma;
    double omegaQ;
    double eta;
};

class cSPINDATA
{
public:
    cSPINDATA();
    SpinProperty getData(string name) {return data[name]; };
private:
    map<string, SpinProperty> data;

};
#endif
