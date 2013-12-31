/*
 * =====================================================================================
 *
 *       Filename:  matrix.cu
 *
 *    Description:  Use the most optimized log&exp method.
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

const int gf_width = 8;
const int field_size = 1 << gf_width;

__device__ __const__ uint16_t gfexp[1021] = { 1,  2,  4,  8,  16,  32,  64,  128,  29,  58,  116,  232,  205,  135,  19,  38,  76,  152,  45,  90,  180,  117,  234,  201,  143,  3,  6,  12,  24,  48,  96,  192,  157,  39,  78,  156,  37,  74,  148,  53,  106,  212,  181,  119,  238,  193,  159,  35,  70,  140,  5,  10,  20,  40,  80,  160,  93,  186,  105,  210,  185,  111,  222,  161,  95,  190,  97,  194,  153,  47,  94,  188,  101,  202,  137,  15,  30,  60,  120,  240,  253,  231,  211,  187,  107,  214,  177,  127,  254,  225,  223,  163,  91,  182,  113,  226,  217,  175,  67,  134,  17,  34,  68,  136,  13,  26,  52,  104,  208,  189,  103,  206,  129,  31,  62,  124,  248,  237,  199,  147,  59,  118,  236,  197,  151,  51,  102,  204,  133,  23,  46,  92,  184,  109,  218,  169,  79,  158,  33,  66,  132,  21,  42,  84,  168,  77,  154,  41,  82,  164,  85,  170,  73,  146,  57,  114,  228,  213,  183,  115,  230,  209,  191,  99,  198,  145,  63,  126,  252,  229,  215,  179,  123,  246,  241,  255,  227,  219,  171,  75,  150,  49,  98,  196,  149,  55,  110,  220,  165,  87,  174,  65,  130,  25,  50,  100,  200,  141,  7,  14,  28,  56,  112,  224,  221,  167,  83,  166,  81,  162,  89,  178,  121,  242,  249,  239,  195,  155,  43,  86,  172,  69,  138,  9,  18,  36,  72,  144,  61,  122,  244,  245,  247,  243,  251,  235,  203,  139,  11,  22,  44,  88,  176,  125,  250,  233,  207,  131,  27,  54,  108,  216,  173,  71,  142,  1,  2,  4,  8,  16,  32,  64,  128,  29,  58,  116,  232,  205,  135,  19,  38,  76,  152,  45,  90,  180,  117,  234,  201,  143,  3,  6,  12,  24,  48,  96,  192,  157,  39,  78,  156,  37,  74,  148,  53,  106,  212,  181,  119,  238,  193,  159,  35,  70,  140,  5,  10,  20,  40,  80,  160,  93,  186,  105,  210,  185,  111,  222,  161,  95,  190,  97,  194,  153,  47,  94,  188,  101,  202,  137,  15,  30,  60,  120,  240,  253,  231,  211,  187,  107,  214,  177,  127,  254,  225,  223,  163,  91,  182,  113,  226,  217,  175,  67,  134,  17,  34,  68,  136,  13,  26,  52,  104,  208,  189,  103,  206,  129,  31,  62,  124,  248,  237,  199,  147,  59,  118,  236,  197,  151,  51,  102,  204,  133,  23,  46,  92,  184,  109,  218,  169,  79,  158,  33,  66,  132,  21,  42,  84,  168,  77,  154,  41,  82,  164,  85,  170,  73,  146,  57,  114,  228,  213,  183,  115,  230,  209,  191,  99,  198,  145,  63,  126,  252,  229,  215,  179,  123,  246,  241,  255,  227,  219,  171,  75,  150,  49,  98,  196,  149,  55,  110,  220,  165,  87,  174,  65,  130,  25,  50,  100,  200,  141,  7,  14,  28,  56,  112,  224,  221,  167,  83,  166,  81,  162,  89,  178,  121,  242,  249,  239,  195,  155,  43,  86,  172,  69,  138,  9,  18,  36,  72,  144,  61,  122,  244,  245,  247,  243,  251,  235,  203,  139,  11,  22,  44,  88,  176,  125,  250,  233,  207,  131,  27,  54,  108,  216,  173,  71,  142,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }; 
__device__ __const__ uint16_t gflog[256] = { 510,  0,  1,  25,  2,  50,  26,  198,  3,  223,  51,  238,  27,  104,  199,  75,  4,  100,  224,  14,  52,  141,  239,  129,  28,  193,  105,  248,  200,  8,  76,  113,  5,  138,  101,  47,  225,  36,  15,  33,  53,  147,  142,  218,  240,  18,  130,  69,  29,  181,  194,  125,  106,  39,  249,  185,  201,  154,  9,  120,  77,  228,  114,  166,  6,  191,  139,  98,  102,  221,  48,  253,  226,  152,  37,  179,  16,  145,  34,  136,  54,  208,  148,  206,  143,  150,  219,  189,  241,  210,  19,  92,  131,  56,  70,  64,  30,  66,  182,  163,  195,  72,  126,  110,  107,  58,  40,  84,  250,  133,  186,  61,  202,  94,  155,  159,  10,  21,  121,  43,  78,  212,  229,  172,  115,  243,  167,  87,  7,  112,  192,  247,  140,  128,  99,  13,  103,  74,  222,  237,  49,  197,  254,  24,  227,  165,  153,  119,  38,  184,  180,  124,  17,  68,  146,  217,  35,  32,  137,  46,  55,  63,  209,  91,  149,  188,  207,  205,  144,  135,  151,  178,  220,  252,  190,  97,  242,  86,  211,  171,  20,  42,  93,  158,  132,  60,  57,  83,  71,  109,  65,  162,  31,  45,  67,  216,  183,  123,  164,  118,  196,  23,  73,  236,  127,  12,  111,  246,  108,  161,  59,  82,  41,  157,  85,  170,  251,  96,  134,  177,  187,  204,  62,  90,  203,  89,  95,  176,  156,  169,  160,  81,  11,  245,  22,  235,  122,  117,  44,  215,  79,  174,  213,  233,  230,  231,  173,  232,  116,  214,  244,  234,  168,  80,  88,  175 }; 

__host__ __device__ uint8_t gf_add(uint8_t a, uint8_t b)
{
	return a ^ b;
}

__host__ __device__ uint8_t gf_sub(uint8_t a, uint8_t b)
{
	return gf_add(a, b);
}

__host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b)
{
	int sum_log;
	sum_log = gflog[a] + gflog[b];
	return gfexp[sum_log];
}

__host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
{
	int sum_log;
	sum_log = gflog[a] + gflog[b];
	return gfexp[sum_log];
}

__host__ __device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b)
{
	uint8_t sum_log;
	while(b)
	{
		if(b & 1)
		{
			sum_log ^= a;
		}
		a = (a << 1) ^ (a & 0x80? 0x1d: 0);
		b >>= 1;
	}
	return sum_log;
}

__host__ __device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
{
	uint8_t sum_log;
	while(b)
	{
		if(b & 1)
		{
			sum_log ^= a;
		}
		a = (a << 1) ^ (a & 0x80? 0x1d: 0);
		b >>= 1;
	}
	return sum_log;
}

__host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b)
{
	int diff_log;
	int gf_max_value = field_size - 1;
	diff_log = gflog[a] + gf_max_value - gflog[b];
	return gfexp[diff_log];
}

__host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
{
	int diff_log;
	int gf_max_value = field_size - 1;
	diff_log = gflog[a] + gf_max_value - gflog[b];
	return gfexp[diff_log];
}

__host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power)
{
	int pow_log = (gflog[a] * power) % (field_size - 1);
	return gfexp[pow_log];
}

__host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power, uint8_t *gflog, uint8_t *gfexp)
{
	int pow_log = (gflog[a] * power) % (field_size - 1);
	return gfexp[pow_log];
}

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

	bx = blockIdx.x;
	do {
// Since we have used (TILE_WIDTH_COL, TILE_WIDTH_ROW) as blockDim, these for loops can be optimized out.
//		for(py = ty; py < TILE_WIDTH_ROW; py += blockDim.y)
//		{
//			for(px = tx; px < TILE_WIDTH_COL; px += blockDim.x)
//			{
				py = ty;
				px = tx;
				row = by*TILE_WIDTH_ROW + py;
				col = bx*TILE_WIDTH_COL + px;
				product[py][px] = 0;
				__syncthreads();
			
				if(row < n && col < m)
				{
					for(int i = 0; i < (int)(ceil((float)p / TILE_DEPTH)); i++)
					{
						int bound = min(p, TILE_DEPTH);
						for(int j = tx; j < bound; j += blockDim.x)
						{
							rowVector[py][j] = A[row*p + i*bound + j];
						}
						for(int j = ty; j < bound; j += blockDim.y)
						{		
							colVector[j][px] = B[col + (i*bound + j)*m];
						}
						__syncthreads();
					
						for(int j = 0; j < bound; j++)
						{
							product[py][px] ^= gf_mul(rowVector[py][j], colVector[j][px]);
				//			dist[py][px] = gf_add(dist[py][px], gf_mul(rowVector[py][j], colVector[j][px]));
						}
						__syncthreads();
					}
					C[row*m+col] = product[py][px];
				}
//			}
//		}
		bx += gridDim.x;
		col = bx*TILE_WIDTH_COL + px;
	} while (col < m);
}

// switch rows if the current row is not the pivot row
__global__ void switch_rows(uint8_t *matrix, uint8_t *result, int rowSrc, int rowDes, int size)
{
    int col = threadIdx.y + blockDim.y * blockIdx.y;
    uint8_t oldMatrixItem;
    uint8_t oldResultItem;

    if( col < size )
    {
        oldMatrixItem = matrix[ index(rowSrc, col, size) ];
        matrix[ index(rowSrc, col, size) ] = matrix[ index(rowDes, col, size) ];
        matrix[ index(rowDes, col, size) ] = oldMatrixItem; 

        oldResultItem = result[ index(rowSrc, col, size) ];
        result[ index(rowSrc, col, size) ] = result[ index(rowDes, col, size) ];
        result[ index(rowDes, col, size) ] = oldResultItem; 
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
        oldMatrixItem = matrix[ index(row, colSrc, size) ];
        matrix[ index(row, colSrc, size) ] = matrix[ index(row, colDes, size) ];
        matrix[ index(row, colDes, size) ] = oldMatrixItem; 

        oldResultItem = result[ index(row, colSrc, size) ];
        result[ index(row, colSrc, size) ] = result[ index(row, colDes, size) ];
        result[ index(row, colSrc, size) ] = oldResultItem; 
    }
} 

// normalize the row by the pivot value
__global__ void normalize_pivot_row(uint8_t *matrix, uint8_t *result, int row, int size)
{
    int ty = threadIdx.y;
	int col = blockDim.y*blockIdx.y + ty;

    __shared__ uint8_t pivotValue;

    if( col < size )
    {
    	// let the first thread of loads the pivotValue
        if ( ty == 0 )
		{
            pivotValue = matrix[ index(row, row, size) ];
		}
        __syncthreads();
		// Normalize the pivot row!
		// Every thread divides the element of its position with the pivotValue
        matrix[ index(row, col, size)] = gf_div(matrix[ index(row, col, size) ], pivotValue);
        result[ index(row, col, size)] = gf_div(result[ index(row, col, size) ], pivotValue);
    }
}

// normalize the column by the pivot value
__global__ void normalize_pivot_col(uint8_t *matrix, uint8_t *result, int col, int size)
{
    int ty = threadIdx.y;
	int row = blockDim.y*blockIdx.y + ty;

    __shared__ uint8_t pivotValue;

    if( col < size )
    {
    	// let the first thread of loads the pivotValue
        if ( ty == 0 )
		{
            pivotValue = matrix[ index(col, col, size) ];
		}
        __syncthreads();
		// Normalize the pivot row!
		// Every thread divides the element of its position with the pivotValue
        matrix[ index(row, col, size)] = gf_div(matrix[ index(row, col, size) ], pivotValue);
        result[ index(row, col, size)] = gf_div(result[ index(row, col, size) ], pivotValue);
    }
}

// eliminate by row to make the pivot column become reduced echelon form
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

    if ( row < size )
    {
        if ( ty == 0 )
        {
            matrixPivotValue = matrix[ index(pivotIndex, col, size) ];
            resultPivotValue = result[ index(pivotIndex, col, size) ];
        }
        pivotCol[ty] = matrix[ index(row, pivotIndex, size) ];
        
        matrixCol[ty] = matrix[ index(row, col, size) ]; 
        resultCol[ty] = result[ index(row, col, size) ]; 
        __syncthreads();

		// substraction in GF
		// make the pivotCol become reduced echelon form
        if ( row != pivotIndex )
        {
			matrix[ index(row, col, size) ] = matrixCol[ty] ^ gf_mul(pivotCol[ty], matrixPivotValue);
			result[ index(row, col, size) ] = resultCol[ty] ^ gf_mul(pivotCol[ty], resultPivotValue);
        }
    }
}

// eliminate by column to make the pivot row become reduced echelon form
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

    if ( row < size )
    {
        if ( ty == 0 )
        {
            matrixPivotValue = matrix[ index(row, pivotIndex, size) ];
            resultPivotValue = result[ index(row, pivotIndex, size) ];
        }
        pivotRow[ty] = matrix[ index(pivotIndex, col, size) ];
        
        matrixCol[ty] = matrix[ index(row, col, size) ]; 
        resultCol[ty] = result[ index(row, col, size) ]; 
        __syncthreads();

		// substraction in GF
		// make the pivotRow become reduced echelon form
        if ( col != pivotIndex )
        {
			matrix[ index(row, col, size) ] = matrixCol[ty] ^ gf_mul(pivotRow[ty], matrixPivotValue);
			result[ index(row, col, size) ] = resultCol[ty] ^ gf_mul(pivotRow[ty], resultPivotValue);
        }
    }
}

// generate an identity matrix
__global__ void get_identity_matrix(uint8_t *result, int size)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.y + threadIdx.y;

    if ( row == col )
	{
        result[ index(row, col, size) ] = 1;
	}
    else
	{
        result[ index(row, col, size) ] = 0;
	}
}

// find the pivot index in the given row/column
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
// Gaussian/Gauss-Jordan elimination
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
	encodingMatrix[i*col + j] = gf_pow(j+1, i);
}

__global__ void encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	matrix_mul(parityCoeff, dataChunk, codeChunk, parityBlockNum, nativeBlockNum, chunkSize);
}

__global__ void decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	matrix_mul(parityCoeff, codeChunk, dataChunk, nativeBlockNum, nativeBlockNum, chunkSize);
}

