#include "include/spin/SpinGrouping.h"

cSpinGrouping::cSpinGrouping(VECT2 connection_matrix)
{
    _connection_matrix=connection_matrix;
    cout << "cSpinGrouping constructor is called." << endl;
}

cSpinGrouing::~cSpinGrouping()
{
    cout << "cSpinGrouping destructor is called." << endl;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinDepthFirstPathTracing
cSpinDepthFirstPathTracing::cSpinDepthFirstPathTracing()
{

    cout << "need a connection matrix." << endl;
}

cSpinDepthFirstPathTracing::cSpinDepthFirstPathTracing(VECT2 connection_matrix)
{
    cout << "cSpinDepthFirstPathTracing constructor is called." << endl;
}

cSpinDepthFirstPathTracing::~cSpinDepthFirstPathTracing()
{
    cout << "cSpinDepthFirstPathTracing destructor is called." << endl;
}

cSpinDepthFirstPathTracing::generate()
{
    
}
