/*
 * =====================================================================================
 *
 *       Filename:  matrix.cu
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/21/2012 07:38:17 PM
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Shuai YUAN (yszheda AT gmail.com), 
 *        Company:  
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "matrix.h"
// #include "gf16.h"

#define DEBUG
// __shared__ uint8_t gflog[256];
// __shared__ uint8_t gfexp[256];

__constant__ uint8_t gfexp[16] = { 1,  2,  4,  8,  3,  6,  12,  11,  5,  10,  7,  14,  15,  13,  9,  0 }; 
__constant__ uint8_t gflog[16] = { 0,  0,  1,  4,  2,  8,  5,  10,  3,  14,  9,  7,  6,  13,  11,  12 }; 
__constant__ uint8_t gfmul_16[16][16] = { { 0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 }, \
 { 0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  9 ,  10 ,  11 ,  12 ,  13 ,  14 ,  15 }, \
 { 0 ,  2 ,  4 ,  6 ,  8 ,  10 ,  12 ,  14 ,  3 ,  1 ,  7 ,  5 ,  11 ,  9 ,  15 ,  13 }, \
 { 0 ,  3 ,  6 ,  5 ,  12 ,  15 ,  10 ,  9 ,  11 ,  8 ,  13 ,  14 ,  7 ,  4 ,  1 ,  2 }, \
 { 0 ,  4 ,  8 ,  12 ,  3 ,  7 ,  11 ,  15 ,  6 ,  2 ,  14 ,  10 ,  5 ,  1 ,  13 ,  9 }, \
 { 0 ,  5 ,  10 ,  15 ,  7 ,  2 ,  13 ,  8 ,  14 ,  11 ,  4 ,  1 ,  9 ,  12 ,  3 ,  6 }, \
 { 0 ,  6 ,  12 ,  10 ,  11 ,  13 ,  7 ,  1 ,  5 ,  3 ,  9 ,  15 ,  14 ,  8 ,  2 ,  4 }, \
 { 0 ,  7 ,  14 ,  9 ,  15 ,  8 ,  1 ,  6 ,  13 ,  10 ,  3 ,  4 ,  2 ,  5 ,  12 ,  11 }, \
 { 0 ,  8 ,  3 ,  11 ,  6 ,  14 ,  5 ,  13 ,  12 ,  4 ,  15 ,  7 ,  10 ,  2 ,  9 ,  1 }, \
 { 0 ,  9 ,  1 ,  8 ,  2 ,  11 ,  3 ,  10 ,  4 ,  13 ,  5 ,  12 ,  6 ,  15 ,  7 ,  14 }, \
 { 0 ,  10 ,  7 ,  13 ,  14 ,  4 ,  9 ,  3 ,  15 ,  5 ,  8 ,  2 ,  1 ,  11 ,  6 ,  12 }, \
 { 0 ,  11 ,  5 ,  14 ,  10 ,  1 ,  15 ,  4 ,  7 ,  12 ,  2 ,  9 ,  13 ,  6 ,  8 ,  3 }, \
 { 0 ,  12 ,  11 ,  7 ,  5 ,  9 ,  14 ,  2 ,  10 ,  6 ,  1 ,  13 ,  15 ,  3 ,  4 ,  8 }, \
 { 0 ,  13 ,  9 ,  4 ,  1 ,  12 ,  8 ,  5 ,  2 ,  15 ,  11 ,  6 ,  3 ,  14 ,  10 ,  7 }, \
 { 0 ,  14 ,  15 ,  1 ,  13 ,  3 ,  2 ,  12 ,  9 ,  7 ,  6 ,  8 ,  4 ,  10 ,  11 ,  5 }, \
 { 0 ,  15 ,  13 ,  2 ,  9 ,  6 ,  4 ,  11 ,  1 ,  14 ,  12 ,  3 ,  8 ,  7 ,  5 ,  10 }, \
}; 

// __host__ __device__ int setup_tables(int w)
// {
// 	unsigned int b;
//    	unsigned int log;
// 	unsigned int x_to_w;
// 	unsigned int prim_poly;
// //	unsigned int r;
// //	unsigned int x;
// //	unsigned int y;
// 
// 	unsigned int prim_poly_4 = 023;
// 	unsigned int prim_poly_8 = 0435;
// 	//uint8_t prim_poly_8 = 285;
// 	unsigned int prim_poly_16 = 0210013;
// 	switch(w) 
// 	{
// 		case 4: prim_poly = prim_poly_4; break;
// 		case 8: prim_poly = prim_poly_8; break;
// 		case 16: prim_poly = prim_poly_16; break;
// 		default: return -1;
// 	}
// 	x_to_w = 1 << w;
// 	b = 1;
// //	r = 0;
// 	for (log = 0; log < x_to_w-1; log++) 
// 	{
// 		/*
// 		r = 0;
// 		x = 1;
// 		y = log;
// 		while(y)
// 		{
// 			printf("y=%d\n",y);
// 			if(y & 1)
// 			{
// 				r = r ^ b;
// 			}
// 			y = y >> 1;
// 			x = x << 1;
// 			if (x & x_to_w) x = x ^ prim_poly;
// 		}
// 			printf("r=%d\n",r);
// 			printf("log=%d\n",log);
// 		*/
// 		if(b>x_to_w) break;
// 		gflog[b] = (uint8_t) log;
// 		gfexp[log] = (uint8_t) b;
// 		b = b << 1;
// 		if (b & x_to_w) 
// 		{
// 			b = b ^ prim_poly;
// 		}
// 	}
// 	return 0;
// }

__host__ __device__ uint8_t gf_add(uint8_t a, uint8_t b)
{
	return a^b;
}

__host__ __device__ uint8_t gf_sub(uint8_t a, uint8_t b)
{
	return gf_add(a, b);
}

__host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b)
{
	uint8_t a_high_nibble = ( a >> 4 ) & 0x0f;
	uint8_t b_high_nibble = ( b >> 4 ) & 0x0f;
	uint8_t result_high_nibble = gfmul_16[a_high_nibble][b_high_nibble];
	uint8_t a_low_nibble = a & 0x0f;
	uint8_t b_low_nibble = b & 0x0f;
	uint8_t result_low_nibble = gfmul_16[a_low_nibble][b_low_nibble];

	return ( (result_high_nibble << 4) | result_low_nibble);
}

// a: a byte composed of two nibbles
// b: a byte that is zero-extended from a nibble
__host__ __device__ uint8_t gf_nibble_mul_byte(uint8_t a, uint8_t b)
{
	uint8_t a_high_nibble = ( a >> 4 ) & 0x0f;
	// b_high_nibble is supposed to be 0
	// uint8_t b_high_nibble = ( b >> 4 ) & 0x0f;
	uint8_t a_low_nibble = a & 0x0f;
	uint8_t b_low_nibble = b & 0x0f;
	uint8_t result_high_nibble = gfmul_16[a_high_nibble][b_low_nibble];
	uint8_t result_low_nibble = gfmul_16[a_low_nibble][b_low_nibble];

	return ( (result_high_nibble << 4) | result_low_nibble);
}

__host__ __device__ uint8_t gf_byte_mul_nibble(uint8_t a, uint8_t b)
{
//	return gf_nibble_mul_byte(b, a);

	// a_high_nibble is supposed to be 0
	// uint8_t a_high_nibble = ( a >> 4 ) & 0x0f;
	uint8_t a_low_nibble = a & 0x0f;
	uint8_t b_low_nibble = b & 0x0f;
	uint8_t b_high_nibble = ( b >> 4 ) & 0x0f;
	uint8_t result_high_nibble = gfmul_16[a_low_nibble][b_high_nibble];
	uint8_t result_low_nibble = gfmul_16[a_low_nibble][b_low_nibble];

	return ( (result_high_nibble << 4) | result_low_nibble);
}

// each of a and b is a byte that is zero-extended from a nibble
__host__ __device__ uint8_t gf_byte_mul_byte(uint8_t a, uint8_t b)
{
	// a_high_nibble is supposed to be 0
	// uint8_t a_high_nibble = ( a >> 4 ) & 0x0f;
	// b_high_nibble is supposed to be 0
	// uint8_t b_high_nibble = ( b >> 4 ) & 0x0f;
	// uint8_t a_low_nibble = a & 0x0f;
	// uint8_t b_low_nibble = b & 0x0f;
	// uint8_t result_high_nibble = gfmul_16[a_high_nibble][b_low_nibble];
	// uint8_t result_low_nibble = gfmul_16[a_low_nibble][b_low_nibble];

	// return (result_high_nibble | result_low_nibble);

	return gfmul_16[a][b];
}

// __host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
// {
// 	int sum_log;
// 	if (a == 0 || b == 0)
// 	{
// 		return 0;
// 	}
// 	sum_log = gflog[a] + gflog[b];
// 	if (sum_log >= NW-1)
// 	{	
// 		sum_log -= NW-1;
// 	}
// 	return gfexp[sum_log];
// }
// 
// __host__ __device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b)
// {
// 	uint8_t sum_log;
// 	while(b)
// 	{
// 		if(b & 1)
// 		{
// 			sum_log ^= a;
// 		}
// 		a = (a << 1) ^ (a & 0x80? 0x1d: 0);
// 		b >>= 1;
// 	}
// 	return sum_log;
// }
// 
// __host__ __device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
// {
// 	uint8_t sum_log;
// 	while(b)
// 	{
// 		if(b & 1)
// 		{
// 			sum_log ^= a;
// 		}
// 		a = (a << 1) ^ (a & 0x80? 0x1d: 0);
// 		b >>= 1;
// 	}
// 	return sum_log;
// }

// __device__ uint8_t gf_div_nibble(uint8_t a_nibble, uint8_t b_nibble)
inline void gf_div_nibble(uint8_t a_nibble, uint8_t b_nibble, int *result_nibble)
{
	int diff_log;
	if (a_nibble == 0)
	{	
//		return 0;
		*(result_nibble) = 0;
	}
	/* Can’t divide by 0 */
	if (b_nibble == 0)
	{
//		return -1;
		*(result_nibble) = -1;
	}
	diff_log = gflog[a_nibble] - gflog[b_nibble];
	if (diff_log < 0)
	{	
		diff_log += NW-1;
	}
//	return gfexp[diff_log];
	*(result_nibble) = gfexp[diff_log];
}

__host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b)
{
	uint8_t a_high_nibble = ( a >> 4 ) & 0x0f;
	uint8_t b_high_nibble = ( b >> 4 ) & 0x0f;

	uint8_t a_low_nibble = a & 0x0f;
	uint8_t b_low_nibble = b & 0x0f;

//	int result_high_nibble = gf_div_nibble(a_high_nibble, b_high_nibble);
//	int result_low_nibble = gf_div_nibble(a_low_nibble, b_low_nibble);

	int result_high_nibble;
	int result_low_nibble;
	gf_div_nibble(a_high_nibble, b_high_nibble, &result_high_nibble);
	gf_div_nibble(a_low_nibble, b_low_nibble, &result_low_nibble);

	/* Can’t divide by 0 */
	if (result_high_nibble == -1 | result_low_nibble == -1)
	{
		return -1;
	}

	return ( (result_high_nibble << 4) | result_low_nibble);
}

// each of a and b is a byte that is zero-extended from a nibble
__host__ __device__ uint8_t gf_byte_div_byte(uint8_t a, uint8_t b)
{
	int diff_log;
	if (a == 0)
	{	
		return 0;
	}
	/* Can’t divide by 0 */
	if (b == 0)
	{
		return -1;
	}
	diff_log = gflog[a] - gflog[b];
	if (diff_log < 0)
	{	
		diff_log += NW-1;
	}

	assert( gfmul_16[ gfexp[diff_log] ][b] == a );
		
	return gfexp[diff_log];
	
}

// __host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
// {
// 	int diff_log;
// 	if (a == 0)
// 	{	
// 		return 0;
// 	}
// 	/* Can’t divide by 0 */
// 	if (b == 0)
// 	{
// 		return -1;
// 	}
// //	diff_log = (gflog[a] - gflog[b]) % (NW-1);
// 	diff_log = gflog[a] - gflog[b];
// 	if (diff_log < 0)
// 	{	
// 		diff_log += NW-1;
// 	}
// 	return gfexp[diff_log];
// }

__host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power)
{
	uint8_t a_high_nibble = ( a >> 4 ) & 0x0f;
	uint8_t a_low_nibble = a & 0x0f;
	int pow_log_high_nibble = (gflog[a_high_nibble] * power) % (NW-1);
	int pow_log_low_nibble = (gflog[a_low_nibble] * power) % (NW-1);
	return ( (gfexp[pow_log_high_nibble] << 4) | gfexp[pow_log_low_nibble]);
}

// each of a and power is a byte that is zero-extended from a nibble
__host__ __device__ uint8_t gf_byte_pow_byte(uint8_t a, uint8_t power)
{
	int pow_log = (gflog[a] * power) % (NW-1);
	return gfexp[pow_log];
}

// __host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power, uint8_t *gflog, uint8_t *gfexp)
// {
// 	uint8_t a_high_nibble = ( a >> 4 ) & 0x0f;
// 	uint8_t a_low_nibble = a & 0x0f;
// 	int pow_log_high_nibble = (gflog[a_high_nibble] * power) % (NW-1);
// 	int pow_log_low_nibble = (gflog[a_low_nibble] * power) % (NW-1);
// 	return (gfexp[pow_log_high_nibble] | gfexp[pow_log_low_nibble]);
// }

// input matrix A and B, compute the product matrix C=AB
// A: nxp
// B: pxm
// C: nxm
__device__ void matrix_mul(unsigned char *A, unsigned char *B, unsigned char *C, int n, int p, int m)
{
	__shared__ int rowVector[TILE_WIDTH_ROW][TILE_DEPTH];
	__shared__ int colVector[TILE_DEPTH][TILE_WIDTH_COL];
	__shared__ int product[TILE_WIDTH_ROW][TILE_WIDTH_COL];

	int bx = blockIdx.x;
   	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int row;
	int col;
	int px;
	int py;	

//	setup_tables(8);
	__syncthreads();

	for(bx=blockIdx.x; bx< (int)(ceil((float)m/gridDim.x)); bx+=gridDim.x )
	{
		for(py=ty; py<TILE_WIDTH_ROW; py+=blockDim.y)
		{
			for(px=tx; px<TILE_WIDTH_COL; px+=blockDim.x)
			{
				row = by*TILE_WIDTH_ROW+py;
				col = bx*TILE_WIDTH_COL+px;
				product[py][px] = 0;
				__syncthreads();
				
				if(row < n && col < m)
				{
					for(int i=0; i<(int)(ceil((float)p/TILE_DEPTH)); i++)
					{
						int bound = min(p, TILE_DEPTH);
		/*
						for(int j=tx; j<TILE_DEPTH; j+=blockDim.x)
						{
							rowVector[py][j] = A[row*p+i*TILE_DEPTH+j];
						}
						for(int j=ty; j<TILE_DEPTH; j+=blockDim.y)
						{		
							colVector[j][px] = B[col+(i*TILE_DEPTH+j)*m];
						}
						__syncthreads();
					
						for(int j=0; j<TILE_DEPTH; j++)
						{
							product[py][px] ^= gf_mul(rowVector[py][j], colVector[j][px]);
		//					dist[py][px] = gf_add(dist[py][px], gf_mul(rowVector[py][j], colVector[j][px]));
						}
		*/
						for(int j=tx; j<bound; j+=blockDim.x)
						{
							rowVector[py][j] = A[row*p+i*bound+j];
						}
						for(int j=ty; j<bound; j+=blockDim.y)
						{		
							colVector[j][px] = B[col+(i*bound+j)*m];
						}
						__syncthreads();
					
						for(int j=0; j<bound; j++)
						{
//							product[py][px] ^= gf_mul(rowVector[py][j], colVector[j][px]);
							product[py][px] ^= gf_byte_mul_nibble(rowVector[py][j], colVector[j][px]);

		//					dist[py][px] = gf_add(dist[py][px], gf_mul(rowVector[py][j], colVector[j][px]));
						}
						__syncthreads();
					}
					C[row*m+col] = product[py][px];
				}
			}
		}
	}
}

// switch rows if the current row is not the pivot row
__global__ void switch_rows(uint8_t *matrix, uint8_t *result, int rowSrc, int rowDes, int size)
{
    int col = threadIdx.y + blockDim.y * blockIdx.y;
    uint8_t oldMatrixItem;
    uint8_t oldResultItem;

    if( col < size )
    {
        oldMatrixItem = matrix[ IDC2D(rowSrc, col, size) ];
        matrix[ IDC2D(rowSrc, col, size) ] = matrix[ IDC2D(rowDes, col, size) ];
        matrix[ IDC2D(rowDes, col, size) ] = oldMatrixItem; 

        oldResultItem = result[ IDC2D(rowSrc, col, size) ];
        result[ IDC2D(rowSrc, col, size) ] = result[ IDC2D(rowDes, col, size) ];
        result[ IDC2D(rowDes, col, size) ] = oldResultItem; 
    }
} 
// switch columns if the current row is not the pivot row
__global__ void switch_columns(uint8_t *matrix, uint8_t *result, int colSrc, int colDes, int size)
{
    int row = threadIdx.y + blockDim.y * blockIdx.y;
    uint8_t oldMatrixItem;
    uint8_t oldResultItem;

    if( row < size )
    {
        oldMatrixItem = matrix[ IDC2D(row, colSrc, size) ];
        matrix[ IDC2D(row, colSrc, size) ] = matrix[ IDC2D(row, colDes, size) ];
        matrix[ IDC2D(row, colDes, size) ] = oldMatrixItem; 

        oldResultItem = result[ IDC2D(row, colSrc, size) ];
        result[ IDC2D(row, colSrc, size) ] = result[ IDC2D(row, colDes, size) ];
        result[ IDC2D(row, colSrc, size) ] = oldResultItem; 
    }
} 

// normalize the row by the pivot value
__global__ void normalize_pivot_row(uint8_t *matrix, uint8_t *result, int row, int size)
{
    int ty = threadIdx.y;
	int col = blockDim.y*blockIdx.y + ty;

    __shared__ uint8_t pivotValue;

//	setup_tables(8);
	__syncthreads();

    if( col < size )
    {
    	// let the first thread of loads the pivotValue
        if ( ty == 0 )
		{
            pivotValue = matrix[ IDC2D(row, row, size) ];
		}
        __syncthreads();
	// Normalize the pivot row!
	// Every thread divides the element of its position with the pivotValue
//        matrix[ IDC2D(row, col, size)] = gf_div(matrix[ IDC2D(row, col, size) ], pivotValue);
//        result[ IDC2D(row, col, size)] = gf_div(result[ IDC2D(row, col, size) ], pivotValue);
        matrix[ IDC2D(row, col, size)] = gf_byte_div_byte(matrix[ IDC2D(row, col, size) ], pivotValue);
        result[ IDC2D(row, col, size)] = gf_byte_div_byte(result[ IDC2D(row, col, size) ], pivotValue);
    }
}
// normalize the column by the pivot value
__global__ void normalize_pivot_col(uint8_t *matrix, uint8_t *result, int col, int size)
{
    int ty = threadIdx.y;
	int row = blockDim.y*blockIdx.y + ty;

    __shared__ uint8_t pivotValue;

//	setup_tables(8);
	__syncthreads();

    if( col < size )
    {
    	// let the first thread of loads the pivotValue
        if ( ty == 0 )
		{
            pivotValue = matrix[ IDC2D(col, col, size) ];
		}
        __syncthreads();
	// Normalize the pivot row!
	// Every thread divides the element of its position with the pivotValue
//        matrix[ IDC2D(row, col, size)] = gf_div(matrix[ IDC2D(row, col, size) ], pivotValue);
//        result[ IDC2D(row, col, size)] = gf_div(result[ IDC2D(row, col, size) ], pivotValue);
        matrix[ IDC2D(row, col, size)] = gf_byte_div_byte(matrix[ IDC2D(row, col, size) ], pivotValue);
        result[ IDC2D(row, col, size)] = gf_byte_div_byte(result[ IDC2D(row, col, size) ], pivotValue);
    }
}

//eliminate by row to make the pivot column become reduced echelon form
__global__ void eliminate_by_row(uint8_t *matrix, uint8_t *result, int pivotIndex, int size)
{
    int ty = threadIdx.y;

	int row = blockDim.y * blockIdx.y + threadIdx.y;
	int col = blockIdx.x;

    __shared__ uint8_t pivotCol[ SINGLE_BLOCK_SIZE ];

    __shared__ uint8_t matrixPivotValue;
    __shared__ uint8_t resultPivotValue;
    __shared__ uint8_t matrixCol[ SINGLE_BLOCK_SIZE ];
    __shared__ uint8_t resultCol[ SINGLE_BLOCK_SIZE];

//	setup_tables(8);
	__syncthreads();

    if ( row < size )
    {
        if ( ty == 0 )
        {
            matrixPivotValue = matrix[ IDC2D(pivotIndex, col, size) ];
            resultPivotValue = result[ IDC2D(pivotIndex, col, size) ];
        }
        pivotCol[ty] = matrix[ IDC2D(row, pivotIndex, size) ];
        
        matrixCol[ty] = matrix[ IDC2D(row, col, size) ]; 
        resultCol[ty] = result[ IDC2D(row, col, size) ]; 
        __syncthreads();

		// substraction in GF
		// make the pivotCol become reduced echelon form
        if ( row != pivotIndex )
        {
//			matrix[ IDC2D(row, col, size) ] = matrixCol[ty] ^ gf_mul(pivotCol[ty], matrixPivotValue);
//			result[ IDC2D(row, col, size) ] = resultCol[ty] ^ gf_mul(pivotCol[ty], resultPivotValue);
			matrix[ IDC2D(row, col, size) ] = matrixCol[ty] ^ gf_byte_mul_byte(pivotCol[ty], matrixPivotValue);
			result[ IDC2D(row, col, size) ] = resultCol[ty] ^ gf_byte_mul_byte(pivotCol[ty], resultPivotValue);
        }
    }
}

//eliminate by column to make the pivot row become reduced echelon form
__global__ void eliminate_by_col(uint8_t *matrix, uint8_t *result, int pivotIndex, int size)
{
    int ty = threadIdx.y;

	int row = blockIdx.x;
	int col = blockDim.y * blockIdx.y + threadIdx.y;

    __shared__ uint8_t pivotRow[ SINGLE_BLOCK_SIZE ];

    __shared__ uint8_t matrixPivotValue;
    __shared__ uint8_t resultPivotValue;
    __shared__ uint8_t matrixCol[ SINGLE_BLOCK_SIZE ];
    __shared__ uint8_t resultCol[ SINGLE_BLOCK_SIZE];

//	setup_tables(8);
	__syncthreads();

    if ( row < size )
    {
        if ( ty == 0 )
        {
            matrixPivotValue = matrix[ IDC2D(row, pivotIndex, size) ];
            resultPivotValue = result[ IDC2D(row, pivotIndex, size) ];
        }
        pivotRow[ty] = matrix[ IDC2D(pivotIndex, col, size) ];
        
        matrixCol[ty] = matrix[ IDC2D(row, col, size) ]; 
        resultCol[ty] = result[ IDC2D(row, col, size) ]; 
        __syncthreads();

		// substraction in GF
		// make the pivotRow become reduced echelon form
        if ( col != pivotIndex )
        {
//			matrix[ IDC2D(row, col, size) ] = matrixCol[ty] ^ gf_mul(pivotRow[ty], matrixPivotValue);
//			result[ IDC2D(row, col, size) ] = resultCol[ty] ^ gf_mul(pivotRow[ty], resultPivotValue);
			matrix[ IDC2D(row, col, size) ] = matrixCol[ty] ^ gf_byte_mul_byte(pivotRow[ty], matrixPivotValue);
			result[ IDC2D(row, col, size) ] = resultCol[ty] ^ gf_byte_mul_byte(pivotRow[ty], resultPivotValue);
        }
    }
}

//generate an identity matrix
__global__ void get_identity_matrix(uint8_t *result, int size)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.y + threadIdx.y;

    if ( row == col )
	{
        result[ IDC2D(row, col, size) ] = 1;
	}
    else
	{
        result[ IDC2D(row, col, size) ] = 0;
	}
}

//find the pivot index in the given row/column
int get_pivot_index(uint8_t *vector, int index, int size)
{
    int pivotIndex = -1;
    int i = index;
    while( pivotIndex == -1 && i < size )
    {
        pivotIndex = (vector[i] > 0)? i: -1;        
        i++;
    }
    return pivotIndex;
}

#ifdef DEBUG
void show_squre_matrix_debug(uint8_t *matrix, int size)
{
	int i;
	int j;
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			printf("%d ", matrix[i*size+j]);
		}
		printf("\n");
	}
		printf("\n");
}
#endif

// compute the inverse of a given matrix
// Gaussian/Gauss–Jordan elimination
extern "C"
void invert_matrix(uint8_t *matrix_dev, uint8_t *result_dev, int size)
{
	int row;
	int pivotIndex;
    uint8_t currentRow[size];
    int currentRowSize = size*sizeof(uint8_t);

    dim3 gimGrid( (int)(ceil( (float)size / SQUARE_BLOCK_SIZE)), (int)(ceil( (float)size / SQUARE_BLOCK_SIZE)) );
    dim3 gimBlock( min(size, SQUARE_BLOCK_SIZE), min(size, SQUARE_BLOCK_SIZE) );
    get_identity_matrix<<< gimGrid, gimBlock >>>(result_dev, size);
//	cudaDeviceSynchronize();
	
	for( row = 0; row < size; row++ )
    {
		// check whether the leading coefficient of the current row is in the 'index'th column
		int index = row;
        cudaMemcpy(currentRow, matrix_dev+row*size, currentRowSize, cudaMemcpyDeviceToHost);
        pivotIndex = get_pivot_index(currentRow, index, size);
        if( pivotIndex != row )
		{
			dim3 scGrid(1, (int)(ceil( (float)size / SINGLE_BLOCK_SIZE )));
			dim3 scBlock(1, min(size, SINGLE_BLOCK_SIZE)); 
            switch_columns<<< scGrid, scBlock >>>(matrix_dev, result_dev, index, pivotIndex, size);
		}
		cudaDeviceSynchronize();

		dim3 nprGrid(1, (int)(ceil( (float)size / SINGLE_BLOCK_SIZE )));
		dim3 nprBlock(1, min(size, SINGLE_BLOCK_SIZE)); 
    	// Normalize the pivot row
        normalize_pivot_row<<< nprGrid, nprBlock >>>(matrix_dev, result_dev, index, size);
//    	// Normalize the pivot column
//        normalize_pivot_col<<< nprGrid, linearBlock >>>(matrix_dev, result_dev, index, size);
		cudaDeviceSynchronize();

		dim3 ebrGrid(size, (int)(ceil( (float)size / SINGLE_BLOCK_SIZE )));
		dim3 ebrBlock(1, min(size, SINGLE_BLOCK_SIZE)); 
        eliminate_by_row<<< ebrGrid, ebrBlock >>>(matrix_dev, result_dev, row, size);
		cudaDeviceSynchronize();

#ifdef DEBUG
uint8_t matrix_host[size*size];
cudaMemcpy(matrix_host, matrix_dev, size*size, cudaMemcpyDeviceToHost);
printf("matrix:\n");
show_squre_matrix_debug(matrix_host, size);
uint8_t result_host[size*size];
cudaMemcpy(result_host, result_dev, size*size, cudaMemcpyDeviceToHost);
printf("result:\n");
show_squre_matrix_debug(result_host, size);
#endif
    }

}

__global__ void gen_encoding_matrix(uint8_t *encodingMatrix, int row, int col)
{
	int i = threadIdx.x;
	int j = threadIdx.y;
//	setup_tables(8);
	__syncthreads();
//	encodingMatrix[i*col + j] = gf_pow(j+1, i);
	encodingMatrix[i*col + j] = gf_byte_pow_byte(j+1, i);
}

__global__ void encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	matrix_mul(parityCoeff, dataChunk, codeChunk, parityBlockNum, nativeBlockNum, chunkSize);
/*
	int currentSize = chunkSize;
	for(int i=0; i<(int)(ceil((float)chunkSize/SINGLE_GRID_SIZE)); i++)
	{
		if(chunkSize-(i+1)*SINGLE_GRID_SIZE < 0)
		{
			currentSize = chunkSize - i*SINGLE_GRID_SIZE;
		}
		matrix_mul(parityCoeff, dataChunk+i*SINGLE_GRID_SIZE, codeChunk+i*SINGLE_GRID_SIZE, parityBlockNum, nativeBlockNum, currentSize);
	}
*/
}

__global__ void decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	matrix_mul(parityCoeff, codeChunk, dataChunk, nativeBlockNum, nativeBlockNum, chunkSize);
/*
	int currentSize = chunkSize;
	for(int i=0; i<(int)(ceil((float)chunkSize/SINGLE_GRID_SIZE)); i++)
	{
		if(chunkSize-(i+1)*SINGLE_GRID_SIZE < 0)
		{
			currentSize = chunkSize - i*SINGLE_GRID_SIZE;
		}
		matrix_mul(parityCoeff, codeChunk+i*SINGLE_GRID_SIZE, dataChunk+i*SINGLE_GRID_SIZE, parityBlockNum, nativeBlockNum, currentSize);
	}
*/
}




