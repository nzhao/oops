#ifndef XMLREADER_H
#define XMLREADER_H
#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "include/rapidxml-1.13/rapidxml.hpp"
#include <iomanip>
#include <armadillo>

using namespace rapidxml;
using namespace std;
using namespace arma;

typedef map<pair<string, string>, pair<string, string> > PARA_MAP;
const int     DEFAULT_INT_PARAMETER = 0;
const double  DEFAULT_DOUBLE_PARAMETER = 0.0;
const string  DEFAULT_STRING_PARAMETER = "None";

class ConfigXML
{
public:
    ConfigXML() {};
    ConfigXML(const ConfigXML& cfg) {_parameters = cfg._parameters;};
    ConfigXML(string filename);
    ~ConfigXML() {};

    void   printParameters() const ;
    int    getIntParameter(string section_name, string para_name) const;
    double getDoubleParameter(string section_name, string para_name) const;
    string getStringParameter(string section_name, string para_name) const;
    vec    getVectorParameter(string section_name, string para_name) const;
    PARA_MAP getParameters() const {return _parameters;};

protected:
private:
    mutable PARA_MAP       _parameters;
};
#endif
