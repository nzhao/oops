#ifndef MISC_H
#define MISC_H
#include <armadillo>

template<class T> 
double distance(T& obj1, T& obj2) { return norm(obj1.get_coordinate() - obj2.get_coordinate()); };
#endif
