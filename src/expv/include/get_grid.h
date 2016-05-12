//#include "cuda_runtime.h"

int get_block( const unsigned long nthreads, const dim3 block_limit, dim3 *block_size );
int get_grid( const unsigned long nthreads, const unsigned int nthreads_max_per_block, const dim3 grid_limit, const dim3 block_limit, dim3 *grid_size, dim3 *block_size );
