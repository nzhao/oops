#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

typedef vector<int> INDICES;
typedef vector< cx_mat > TERM;
typedef double MULTIPLIER;

typedef vector< INDICES > INDEX_LIST;
typedef vector< vector<TERM> > MAT_LIST;
typedef vector< Col<MULTIPLIER>  > COEFF_LIST;

typedef vector<int> DIM_LIST;

