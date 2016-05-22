#include "include/misc/print_program_options.h"

using namespace std;
namespace po = boost::program_options;

void PrintVariableMap(po::variables_map vm) 
{
    cout << "#############################################################" << endl;
    cout << "Parameters Using Default Values: " << endl;
    cout << "-------------------------------------------------------------" << endl;
    for (po::variables_map::iterator it = vm.begin(); it != vm.end(); ++it) 
        if (vm[it->first].defaulted() || it->second.defaulted()) 
        {
            int isString, isInt, isDouble;
            cout << std::right << std::setw(15) << it->first << "  =  " << setw(30) << std::left;
            try { cout << vm[it->first].as<string>() << " (string)" << std::endl; } 
            catch (const boost::bad_any_cast &) { isString = false; }

            try { cout << vm[it->first].as<int>() << " (int)" << std::endl; } 
            catch (const boost::bad_any_cast &) { isInt = false; }

            try { cout << vm[it->first].as<double>() << " (double)" << std::endl; } 
            catch (const boost::bad_any_cast &) { isDouble = false; }
        }

    cout << "-------------------------------------------------------------" << endl;
    cout << "User Defined Parameters: " << endl;
    cout << "-------------------------------------------------------------" << endl;
    for (po::variables_map::iterator it = vm.begin(); it != vm.end(); ++it) 
        if (! vm[it->first].defaulted() && ! it->second.defaulted() ) 
        {
            int isString, isInt, isDouble;
            cout << std::right << setw(15) << it->first << "  =  " << setw(30) << std::left;
            try { cout << vm[it->first].as<string>() << " (string)" << std::endl; } 
            catch (const boost::bad_any_cast &) { isString = false; }

            try { cout << vm[it->first].as<int>() << " (int)" << std::endl; } 
            catch (const boost::bad_any_cast &) { isInt = false; }

            try { cout << vm[it->first].as<double>() << " (double)" << std::endl; } 
            catch (const boost::bad_any_cast &) { isDouble = false; }
        }

    cout << "#############################################################" << endl;
    cout << endl;
}
