//======================================================================
/*
  get_block,  assign threads into a 3D block.
    
    input:
      
      nthreads,     number of threads.
      block_limit,  max dimension of the block.
    
    output:
      
      block_size,   size of the block.
      return value, 0 for all threads are assigned to the block;
                    1 for part of the threads are assigned.
    
  get_grid,   assign threads into a 3D grid, which consists of blocks.
    
    input:
      
      nthreads,     number of threads to assign.
      nthreads_max_per_block, max threads per block.
      grid_limit,   max dimension of the grid.
      block_limit,  max dimension of the block in the grid.
      
    output:
      
      grid_size,    size of the grid.
      block_size,   size of the block.
      return value, 0 for all threads are assigned to the grid;
                    1 for part of the threads are assigned.
*/
//======================================================================

//#include <cstddef>
//#include <cmath>
#include "cuda_runtime.h"
//#include "cublas_v2.h"

//// head;
//int get_block( const unsigned long nthreads, const dim3 block_limit, dim3 *block_size );
//int get_grid( const unsigned long nthreads, const unsigned int nthreads_max_per_block, const dim3 grid_limit, const dim3 block_limit, dim3 *grid_size, dim3 *block_size );

int get_block( const unsigned long nthreads, const dim3 block_limit, dim3 *block_size )
{
  unsigned long   v, v_xy, v_xyz;
  unsigned int    zlimit, ylimit, z, y, x;
  
  int             block_status(1);
  unsigned long   nthreads_max(0);
  dim3            max_size(0,0,0);
  
  v = (unsigned long)block_limit.x * (unsigned long)block_limit.y * (unsigned long)block_limit.z;
  if ( v > nthreads )
    v = nthreads;
  
  zlimit = ( v - 1 ) / ( (unsigned long)block_limit.x * (unsigned long)block_limit.y ) + 1;
  for ( z = ( ( 1 > zlimit ) ? 1 : zlimit ); z <= block_limit.z; z++ )
  {
    v_xy = v / z;
    ylimit = ( v_xy - 1 ) / (unsigned long)block_limit.x + 1;
    for ( y = ( ( 1 > ylimit ) ? 1 : ylimit ); y <= block_limit.y; y++ )
    {
      x = v_xy / y;
      if ( x < 1 ) continue;
      
      v_xyz = (unsigned long)x * (unsigned long)y * (unsigned long)z;
      if ( v_xyz > nthreads_max )
      {
        nthreads_max = v_xyz;
        max_size.x = x;
        max_size.y = y;
        max_size.z = z;
        
        if ( v_xyz == v )
        {
          block_status = 0;
          break;
        }
      }
    }
    if ( block_status == 0 )
      break;
  }
  
  if ( nthreads_max == nthreads )
    block_status = 0;
  else
    block_status = 1;
  (*block_size) = max_size;
  
  return block_status;
}

int get_grid( const unsigned long nthreads, const unsigned int nthreads_max_per_block, const dim3 grid_limit, const dim3 block_limit, dim3 *grid_size, dim3 *block_size )
{
  unsigned int    nthreads_per_block_limit, nthreads_per_block;
  unsigned long   nblocks;
  
  int             grid_status(1);
  dim3            bdsize(0,0,0), gdsize(0,0,0);
  unsigned int    nt;
  
  nthreads_per_block_limit = block_limit.x * block_limit.y * block_limit.z;
  if ( nthreads_per_block_limit > nthreads_max_per_block )
    nthreads_per_block_limit = nthreads_max_per_block;
  
  for ( nt = nthreads_per_block_limit; nt >= 1; nt-- )
  {
    get_block( nt, block_limit, &bdsize );
    nthreads_per_block = bdsize.x * bdsize.y * bdsize.z;
    if ( ( nthreads % (unsigned long)nthreads_per_block ) != 0 ) continue;
    
    nblocks = nthreads / (unsigned long)nthreads_per_block;
    grid_status = get_block( nblocks, grid_limit, &gdsize );
    if ( grid_status == 0 ) break;
  }
  
  if ( grid_status == 0 )
  {
    (*grid_size)  = gdsize;
    (*block_size) = bdsize;
  }
  
  return grid_status;
}

//#include <cstddef>
//#include <cmath>
//#include "cuda_runtime.h"
//#include "cublas_v2.h"

//// head;
//int get_block( const size_t *nthreads, const dim3 *block_limit, dim3 *block_size );
//int get_grid( const size_t *nthreads, const size_t *nthreads_max_per_block, const dim3 *grid_limit, const dim3 *block_limit, dim3 *grid_size, dim3 *block_size );

//int get_block( const size_t *nthreads, const dim3 *block_limit, dim3 *block_size )
//{
//  size_t  v, x, y, z, v_xy, v_xyz;
//  size_t  zlimit, ylimit;
//  
//  int     block_status(1);
//  size_t  nthreads_max(0);
//  dim3    max_size(0,0,0);
//  
//  v = (*block_limit).x * (*block_limit).y * (*block_limit).z;
//  if ( v > *nthreads )
//    v = *nthreads;
//  
//  zlimit = ( v - 1 ) / ( (*block_limit).x * (*block_limit).y ) + 1;
//  for ( z = ( ( 1 > zlimit ) ? 1 : zlimit ); z <= (*block_limit).z; z++ )
//  {
//    v_xy = v / z;
//    ylimit = ( v_xy - 1 ) / (*block_limit).x + 1;
//    for ( y = ( ( 1 > ylimit ) ? 1 : ylimit ); y <= (*block_limit).y; y++ )
//    {
//      x = v_xy / y;
//      if ( x < 1 ) continue;
//      
//      v_xyz = x * y * z;
//      if ( v_xyz > nthreads_max )
//      {
//        nthreads_max = v_xyz;
//        max_size.x = x;
//        max_size.y = y;
//        max_size.z = z;
//        
//        if ( v_xyz == v )
//        {
//          block_status = 0;
//          break;
//        }
//      }
//    }
//    if ( block_status == 0 )
//      break;
//  }
//  
//  if ( nthreads_max == *nthreads )
//    block_status = 0;
//  else
//    block_status = 1;
//  (*block_size) = max_size;
//  
//  return block_status;
//}

//int get_grid( const size_t *nthreads, const size_t *nthreads_max_per_block, const dim3 *grid_limit, const dim3 *block_limit, dim3 *grid_size, dim3 *block_size )
//{
//  size_t  nthreads_per_block_limit, nthreads_per_block, nblocks;
//  
//  int     grid_status(1);
//  dim3    bdsize(0,0,0), gdsize(0,0,0);
//  size_t  nt;
//  
//  nthreads_per_block_limit = (*block_limit).x * (*block_limit).y * (*block_limit).z;
//  if ( nthreads_per_block_limit > *nthreads_max_per_block )
//    nthreads_per_block_limit = *nthreads_max_per_block;
//  
//  for ( nt = nthreads_per_block_limit; nt >= 1; nt-- )
//  {
//    get_block( &nt, block_limit, &bdsize );
//    nthreads_per_block = bdsize.x * bdsize.y * bdsize.z;
//    if ( ( *nthreads % nthreads_per_block ) != 0 ) continue;
//    
//    nblocks = *nthreads / nthreads_per_block;
//    grid_status = get_block( &nblocks, grid_limit, &gdsize );
//    if ( grid_status == 0 ) break;
//  }
//  
//  if ( grid_status == 0 )
//  {
//    (*grid_size)  = gdsize;
//    (*block_size) = bdsize;
//  }
//  
//  return grid_status;
//}

