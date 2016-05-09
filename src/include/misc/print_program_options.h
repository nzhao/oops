#include <string.h>
#include <iostream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/any.hpp"

using namespace std;
namespace po = boost::program_options;

typedef std::_Rb_tree_const_iterator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, boost::program_options::variable_value> > bpoIter;
inline void print_iterm(const boost::program_options::variables_map& vm, bpoIter& it)
{/*{{{*/
    cout << std::right << setw(15) << it->first;
    cout << "  =  " << setw(30) << std::left;
    if (((boost::any)it->second.value()).empty()) 
        cout << "(empty)";

    bool is_char;
    try {
        boost::any_cast<const char *>(it->second.value());
        is_char = true;
    } catch (const boost::bad_any_cast &) {
        is_char = false;
    }

    bool is_str;
    try {
        boost::any_cast<std::string>(it->second.value());
        is_str = true;
    } catch (const boost::bad_any_cast &) {
        is_str = false;
    }

    if (((boost::any)it->second.value()).type() == typeid(int)) {
        cout << vm[it->first].as<int>() << " (int)" << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(bool)) {
        cout << vm[it->first].as<bool>() << " (bool)" << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(double)) {
        cout << vm[it->first].as<double>() << " (double)" << std::endl;
    } else if (is_char) {
        cout << vm[it->first].as<const char * >() << std::endl;
    } else if (is_str) {
        std::string temp = vm[it->first].as<std::string>();
        if (temp.size()) {
            cout << temp << " (string)" <<  std::endl;
        } else {
            cout << "true" << std::endl;
        }
    } else { // Assumes that the only remainder is vector<string>
        try {
            std::vector<std::string> vect = vm[it->first].as<std::vector<std::string> >();
            uint i = 0;
            for (std::vector<std::string>::iterator oit=vect.begin();
                    oit != vect.end(); oit++, ++i) {
                cout << "\r> " << it->first << "[" << i << "]=" << (*oit) << std::endl;
            }
        } catch (const boost::bad_any_cast &) {
            cout << "UnknownType(" << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
        }
    }
}/*}}}*/

inline void PrintVariableMap(const boost::program_options::variables_map vm) {
    bpoIter it;
    cout << "#############################################################" << endl;
    for (it = vm.begin(); it != vm.end(); ++it) 
        if (vm[it->first].defaulted() || it->second.defaulted()) 
            print_iterm(vm, it);

    cout << endl;
    for (it = vm.begin(); it != vm.end(); ++it) 
        if (! vm[it->first].defaulted() && ! it->second.defaulted() ) 
            print_iterm(vm, it);
    cout << "#############################################################" << endl;
    cout << endl;
}
