#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

typedef vector<size_t> INDICES; ///< InDICES is a set of spin index 
typedef vector< cx_mat > TERM; ///< a TERM is a list of matrix
typedef double MULTIPLIER; ///< a MULTIPLIER is a number

typedef vector< INDICES > INDEX_LIST; ///< INDEX_LIST is a list of INDICES
typedef vector< vector<TERM> > MAT_LIST; ///< MAT_LIST is ...
typedef vector< Col<MULTIPLIER>  > COEFF_LIST; ///< COEFF_LIST is a ...

typedef vector<size_t> DIM_LIST; ///< DIM_LIST is a list of dim.

