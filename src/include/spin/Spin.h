#ifndef SPIN_H
#define SPIN_H

#include <vector>
#include <string>
#include "include/spin/SpinData.h"

using namespace std;

extern cSPINDATA SPIN_DATABASE;

class cSPIN
{
public:
    cSPIN(vector<double> coord, string isotope_str);
    cSPIN();

    vector<double> get_coordinate() { return coordinate; };
    string get_isotope() { return isotope; };
    int get_multiplicity() { return multiplicity; };
    double get_gamma() { return gamma; };

    void set_coordinate(vector<double> coord) { coordinate = coord; };
    void set_isotope(string iso_str) { isotope =  iso_str; };
    
private:
    vector<double> coordinate;
    string isotope;
    int multiplicity;
    double gamma;
};

#endif
