#ifndef MISC_H
#define MISC_H
#include <armadillo>
#include "include/spin/Spin.h"
#include <cassert>

using namespace arma;

extern cx_double II;

double spin_distance(const cSPIN& spin1, const cSPIN& spin2);

vec dipole(const cSPIN& spin1, const cSPIN& spin2);

vec r_vect(const cSPIN& obj1, const cSPIN& obj2);

vec zeeman(const cSPIN& spin, const vec& magB);

vec dipole_field(const cSPIN& spin, const cSPIN& source_spin, const cx_vec& source_state_vect);

vector<double> Pulse_Timing(string pulsename, int n);

vector<double> Pulse_Interval(string pulsename, int n);

template<class T> vector<T> riffle(T obj1, T obj2, int n)
{
    vector<T> res;

    int q = (n+1)/2;
    for(int i=0; i<q; ++i)
    {
        res.push_back(obj1);
        res.push_back(obj2);
    }

    
    if( n % 2 == 0 )
        res.push_back(obj1);
    return res;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &A, const std::vector<T> &B)
{
    std::vector<T> AB;
    AB.reserve( A.size() + B.size() );                // preallocate memory
    AB.insert( AB.end(), A.begin(), A.end() );        // add A;
    AB.insert( AB.end(), B.begin(), B.end() );        // add B;
    return AB;
}

template <typename T>
std::vector<T> &operator+=(std::vector<T> &A, const std::vector<T> &B)
{
    A.reserve( A.size() + B.size() );                // preallocate memory without erase original data
    A.insert( A.end(), B.begin(), B.end() );         // add B;
    return A;                                        // here A could be named AB
}

template <typename T>
T join_all(const std::vector<T> &list)
{
    T res;
    if( !list.empty() )
    {
        res = list[0];
        for(int i=1; i<list.size(); ++i)
            res += list[i];
    }
    return res;
}
vector<int> base_transform(int num, const vector<int>& base);
int base_number(const vector<int>& num_in_base, const vector<int> base);

template<typename T>
void print_vector(const vector<T>& v)
{
    cout << "[ ";
    for(int i=0; i<v.size(); ++i)
    {
        cout << v[i] ;
        if(i<v.size()-1) 
            cout << ", ";
    }
    cout << " ]" ;
}

#endif
