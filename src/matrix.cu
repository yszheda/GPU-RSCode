/*
 * =====================================================================================
 *
 *       Filename:  matrix.cu
 *
 *    Description:  log&exp method with optimization technique II.
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

const int gfexp_table_size = 1021;
const int gflog_table_size = 256;
__device__ __const__ uint8_t gfexp_cMem[gfexp_table_size] = { 1,  2,  4,  8,  16,  32,  64,  128,  29,  58,  116,  232,  205,  135,  19,  38,  76,  152,  45,  90,  180,  117,  234,  201,  143,  3,  6,  12,  24,  48,  96,  192,  157,  39,  78,  156,  37,  74,  148,  53,  106,  212,  181,  119,  238,  193,  159,  35,  70,  140,  5,  10,  20,  40,  80,  160,  93,  186,  105,  210,  185,  111,  222,  161,  95,  190,  97,  194,  153,  47,  94,  188,  101,  202,  137,  15,  30,  60,  120,  240,  253,  231,  211,  187,  107,  214,  177,  127,  254,  225,  223,  163,  91,  182,  113,  226,  217,  175,  67,  134,  17,  34,  68,  136,  13,  26,  52,  104,  208,  189,  103,  206,  129,  31,  62,  124,  248,  237,  199,  147,  59,  118,  236,  197,  151,  51,  102,  204,  133,  23,  46,  92,  184,  109,  218,  169,  79,  158,  33,  66,  132,  21,  42,  84,  168,  77,  154,  41,  82,  164,  85,  170,  73,  146,  57,  114,  228,  213,  183,  115,  230,  209,  191,  99,  198,  145,  63,  126,  252,  229,  215,  179,  123,  246,  241,  255,  227,  219,  171,  75,  150,  49,  98,  196,  149,  55,  110,  220,  165,  87,  174,  65,  130,  25,  50,  100,  200,  141,  7,  14,  28,  56,  112,  224,  221,  167,  83,  166,  81,  162,  89,  178,  121,  242,  249,  239,  195,  155,  43,  86,  172,  69,  138,  9,  18,  36,  72,  144,  61,  122,  244,  245,  247,  243,  251,  235,  203,  139,  11,  22,  44,  88,  176,  125,  250,  233,  207,  131,  27,  54,  108,  216,  173,  71,  142,  1,  2,  4,  8,  16,  32,  64,  128,  29,  58,  116,  232,  205,  135,  19,  38,  76,  152,  45,  90,  180,  117,  234,  201,  143,  3,  6,  12,  24,  48,  96,  192,  157,  39,  78,  156,  37,  74,  148,  53,  106,  212,  181,  119,  238,  193,  159,  35,  70,  140,  5,  10,  20,  40,  80,  160,  93,  186,  105,  210,  185,  111,  222,  161,  95,  190,  97,  194,  153,  47,  94,  188,  101,  202,  137,  15,  30,  60,  120,  240,  253,  231,  211,  187,  107,  214,  177,  127,  254,  225,  223,  163,  91,  182,  113,  226,  217,  175,  67,  134,  17,  34,  68,  136,  13,  26,  52,  104,  208,  189,  103,  206,  129,  31,  62,  124,  248,  237,  199,  147,  59,  118,  236,  197,  151,  51,  102,  204,  133,  23,  46,  92,  184,  109,  218,  169,  79,  158,  33,  66,  132,  21,  42,  84,  168,  77,  154,  41,  82,  164,  85,  170,  73,  146,  57,  114,  228,  213,  183,  115,  230,  209,  191,  99,  198,  145,  63,  126,  252,  229,  215,  179,  123,  246,  241,  255,  227,  219,  171,  75,  150,  49,  98,  196,  149,  55,  110,  220,  165,  87,  174,  65,  130,  25,  50,  100,  200,  141,  7,  14,  28,  56,  112,  224,  221,  167,  83,  166,  81,  162,  89,  178,  121,  242,  249,  239,  195,  155,  43,  86,  172,  69,  138,  9,  18,  36,  72,  144,  61,  122,  244,  245,  247,  243,  251,  235,  203,  139,  11,  22,  44,  88,  176,  125,  250,  233,  207,  131,  27,  54,  108,  216,  173,  71,  142,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }; 
__device__ __const__ uint16_t gflog_cMem[gflog_table_size] = { 510,  0,  1,  25,  2,  50,  26,  198,  3,  223,  51,  238,  27,  104,  199,  75,  4,  100,  224,  14,  52,  141,  239,  129,  28,  193,  105,  248,  200,  8,  76,  113,  5,  138,  101,  47,  225,  36,  15,  33,  53,  147,  142,  218,  240,  18,  130,  69,  29,  181,  194,  125,  106,  39,  249,  185,  201,  154,  9,  120,  77,  228,  114,  166,  6,  191,  139,  98,  102,  221,  48,  253,  226,  152,  37,  179,  16,  145,  34,  136,  54,  208,  148,  206,  143,  150,  219,  189,  241,  210,  19,  92,  131,  56,  70,  64,  30,  66,  182,  163,  195,  72,  126,  110,  107,  58,  40,  84,  250,  133,  186,  61,  202,  94,  155,  159,  10,  21,  121,  43,  78,  212,  229,  172,  115,  243,  167,  87,  7,  112,  192,  247,  140,  128,  99,  13,  103,  74,  222,  237,  49,  197,  254,  24,  227,  165,  153,  119,  38,  184,  180,  124,  17,  68,  146,  217,  35,  32,  137,  46,  55,  63,  209,  91,  149,  188,  207,  205,  144,  135,  151,  178,  220,  252,  190,  97,  242,  86,  211,  171,  20,  42,  93,  158,  132,  60,  57,  83,  71,  109,  65,  162,  31,  45,  67,  216,  183,  123,  164,  118,  196,  23,  73,  236,  127,  12,  111,  246,  108,  161,  59,  82,  41,  157,  85,  170,  251,  96,  134,  177,  187,  204,  62,  90,  203,  89,  95,  176,  156,  169,  160,  81,  11,  245,  22,  235,  122,  117,  44,  215,  79,  174,  213,  233,  230,  231,  173,  232,  116,  214,  244,  234,  168,  80,  88,  175 }; 
__shared__ uint8_t gfexp[gfexp_table_size];
__shared__ uint16_t gflog[gflog_table_size];

__device__ int setup_tables()
{
//	const int gf_width = 8;
//	const int field_size = 1 << gf_width;
	const unsigned int prim_poly = 0435;
   	int log;
	int exp = 1;
	// use int as book-keeping index instead of unsigned int
	for (log = 0; log < field_size - 1; log++) 
	{
		if(exp > field_size) break;
		gflog[exp] = (uint8_t) log;
		gfexp[log] = (uint8_t) exp;
		if (log < 255)
		{
			gfexp[log + 255] = (uint8_t) exp; 
		}
		exp = exp << 1;
		if (exp & field_size) 
		{
			exp = exp ^ prim_poly;
		}
	}
	int gf_max_value = field_size - 1;
	gflog[0] = 2 * gf_max_value;
	for (log = 2 * gf_max_value; log <= 4 * gf_max_value; ++log)
	{
		gfexp[log] = 0;
	}
	return 0;
}

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
	if (a == 0)
	{	
		return 0;
	}
	// optimize out exception cases
	/*
	// Can't divide by 0
	if (b == 0)
	{
		return -1;
	}
	*/
	const int gf_max_value = (1 << gf_width) - 1;
	diff_log = gflog[a] + gf_max_value - gflog[b];
	return gfexp[diff_log];
}

__host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
{
	int diff_log;
	if (a == 0)
	{	
		return 0;
	}
	// optimize out exception cases
	/*
	// Can't divide by 0
	if (b == 0)
	{
		return -1;
	}
	*/
	const int gf_max_value = (1 << gf_width) - 1;
	diff_log = gflog[a] + gf_max_value - gflog[b];
	return gfexp[diff_log];
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
__global__ void matrix_mul(unsigned char *A, unsigned char *B, unsigned char *C, int n, int p, int m, int tileWidthRow, int tileWidthCol, int tileDepth)
{
	extern __shared__ uint8_t sMem[];
	int rowVectorSize = tileWidthRow * tileDepth;
	int colVectorSize = tileDepth * tileWidthCol;
	int product;

	int bx = blockIdx.x;
   	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int row;
	int col;

//	setup_tables();
//	__syncthreads();

	#pragma unroll
//	for(int j = tx; j < gflog_table_size; j += blockDim.x)
	for(int j = tx * blockDim.y + ty; j < gflog_table_size; j += blockDim.x * blockDim.y)
	{
		gflog[j] = gflog_cMem[j];
	}
	#pragma unroll
//	for(int j = tx; j < gfexp_table_size; j += blockDim.x)
	for(int j = tx * blockDim.y + ty; j < gfexp_table_size; j += blockDim.x * blockDim.y)
	{
		gfexp[j] = gfexp_cMem[j];
	}
	__syncthreads();

	bx = blockIdx.x;
	do {
		row = by * tileWidthRow + ty;
		col = bx * tileWidthCol + tx;
		product = 0;
		__syncthreads();
		
		if(row < n && col < m)
		{
			for(int j = tx; j < tileDepth; j += blockDim.x)
			{
				sMem[ index(ty, j, tileDepth) ] = A[row * p + j];
			}
			for(int j = ty; j < tileDepth; j += blockDim.y)
			{
				sMem[rowVectorSize + index(j, tx, tileWidthCol)] = B[col + j * m];
			}
//			// Since blockDim.x > tileDepth for our applications,
//			// we can fully parallelize loading matrix A into sMem.
//			if (tx < tileDepth)
//			{
//				sMem[ index(ty, tx, tileDepth) ] = A[row * p + tx];
//			}
//			TODO: Assume removing the loop
//			if (ty < tileDepth)
//			{
//				sMem[rowVectorSize + index(ty, tx, tileWidthCol)] = B[col + ty * m];
//			}
			__syncthreads();
			
			for(int j = 0; j < tileDepth; j++)
			{
				product ^= gf_mul(sMem[ index(ty, j, tileDepth) ], sMem[rowVectorSize + index(j, tx, tileWidthCol)]);
			}
			__syncthreads();
			C[row * m + col] = product;
		}
		bx += gridDim.x;
		col = bx * tileWidthCol + tx;
		__syncthreads();
	} while (col < m);
}

// switch rows if the current row is not the pivot row
__global__ void switch_rows(uint8_t *matrix, uint8_t *result, int rowSrc, int rowDes, int size)
{
    int col = threadIdx.y + blockDim.y * blockIdx.y;
    uint8_t oldMatrixItem;
    uint8_t oldResultItem;

    if(col < size)
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

    if(row < size)
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

	setup_tables();
	__syncthreads();

    if(col < size)
    {
    	// let the first thread of loads the pivotValue
        if (ty == 0)
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

	setup_tables();
	__syncthreads();

    if(col < size)
    {
    	// let the first thread of loads the pivotValue
        if (ty == 0)
		{
            pivotValue = matrix[ index(col, col, size) ];
		}
        __syncthreads();
		// Normalize the pivot row!1024
		// Every thread divides the element of its position with the pivotValue
        matrix[ index(row, col, size)] = gf_div(matrix[ index(row, col, size) ], pivotValue);
        result[ index(row, col, size)] = gf_div(result[ index(row, col, size) ], pivotValue);
    }
}

// eliminate by row to make the pivot column become reduced echelon form
__global__ void eliminate_by_row(volatile uint8_t *matrix, volatile uint8_t *result, int pivotIndex, int size)
{
    const int ty = threadIdx.y;
	const int row = blockDim.y * blockIdx.y + threadIdx.y;
	const int col = blockIdx.x;

	setup_tables();
	__syncthreads();

    if (row < size)
    {
		// substraction in GF
		// make the pivotCol become reduced echelon form
        if (row != pivotIndex)
        {
			uint8_t newMatrixValue = matrix[ index(row, col, size) ] ^ gf_mul(matrix[ index(row, pivotIndex, size) ], matrix[ index(pivotIndex, col, size) ]);
			__threadfence();
			__syncthreads();
			uint8_t newResultValue = result[ index(row, col, size) ] ^ gf_mul(matrix[ index(row, pivotIndex, size) ], result[ index(pivotIndex, col, size) ]);
			__threadfence();
			__syncthreads();
			matrix[ index(row, col, size) ] = newMatrixValue;
			__threadfence();
			__syncthreads();
			result[ index(row, col, size) ] = newResultValue;
			__threadfence();
			__syncthreads();
        }
    }
}

// eliminate by column to make the pivot row become reduced echelon form
__global__ void eliminate_by_col(uint8_t *matrix, uint8_t *result, int pivotIndex, int size)
{
    int ty = threadIdx.y;
	int row = blockIdx.x;
	int col = blockDim.y * blockIdx.y + threadIdx.y;

	setup_tables();
	__syncthreads();

    if (col < size)
    {
        if (col != pivotIndex)
        {
			uint8_t newMatrixValue = matrix[ index(row, col, size) ] ^ gf_mul(matrix[ index(pivotIndex, col, size) ], matrix[ index(row, pivotIndex, size) ]);
			__threadfence();
			__syncthreads();
			uint8_t newResultValue = result[ index(row, col, size) ] ^ gf_mul(matrix[ index(pivotIndex, col, size) ], result[ index(row, pivotIndex, size) ]);
			__threadfence();
			__syncthreads();
			matrix[ index(row, col, size) ] = newMatrixValue;
			__threadfence();
			__syncthreads();
			result[ index(row, col, size) ] = newResultValue;
			__threadfence();
			__syncthreads();
        }
    }
}

// generate an identity matrix
__global__ void get_identity_matrix(uint8_t *result, int size)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.y + threadIdx.y;

    if (row == col)
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
    while(pivotIndex == -1 && i < size)
    {
        pivotIndex = (vector[i] > 0)? i: -1;        
        i++;
    }
    return pivotIndex;
}

// #define DEBUG
#ifdef DEBUG
void show_squre_matrix_debug(uint8_t *matrix, int size)
{
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			printf("%d ", matrix[i*size+j]);
		}
		printf("\n");
	}
		printf("\n");
}
#endif

// compute the inverse of a given matrix
// Gaussian elimination
extern "C"
void invert_matrix(uint8_t *matrix_dev, uint8_t *result_dev, int size)
{
	int pivotIndex;
    uint8_t currentRow[size];
    int currentRowSize = size * sizeof(uint8_t);

    dim3 gimGrid((int) (ceil((float) size / SQUARE_BLOCK_SIZE)), (int) (ceil((float) size / SQUARE_BLOCK_SIZE)) );
    dim3 gimBlock(min(size, SQUARE_BLOCK_SIZE), min(size, SQUARE_BLOCK_SIZE));
    get_identity_matrix<<<gimGrid, gimBlock>>>(result_dev, size);
//	cudaDeviceSynchronize();
	
	for(int row = 0; row < size; row++)
    {
		// check whether the leading coefficient of the current row is in the 'index'th column
		int index = row;
        cudaMemcpy(currentRow, matrix_dev + row * size, currentRowSize, cudaMemcpyDeviceToHost);
        pivotIndex = get_pivot_index(currentRow, index, size);
        if(pivotIndex != row)
		{
			dim3 scGrid(1, (int) (ceil((float) size / SINGLE_BLOCK_SIZE)));
			dim3 scBlock(1, min(size, SINGLE_BLOCK_SIZE)); 
            switch_columns<<<scGrid, scBlock>>>(matrix_dev, result_dev, index, pivotIndex, size);
		}
		cudaDeviceSynchronize();

#ifdef DEBUG
uint8_t matrix_host[size * size];
cudaMemcpy(matrix_host, matrix_dev, size * size, cudaMemcpyDeviceToHost);
printf("Current row: %d\n", row);
printf("Step: switch columns\n");
printf("matrix:\n");
show_squre_matrix_debug(matrix_host, size);
uint8_t result_host[size * size];
cudaMemcpy(result_host, result_dev, size * size, cudaMemcpyDeviceToHost);
printf("result:\n");
show_squre_matrix_debug(result_host, size);
#endif

		dim3 nprGrid(1, (int) (ceil((float) size / SINGLE_BLOCK_SIZE)));
		dim3 nprBlock(1, min(size, SINGLE_BLOCK_SIZE)); 
    	// Normalize the pivot row
        normalize_pivot_row<<<nprGrid, nprBlock>>>(matrix_dev, result_dev, index, size);
//    	// Normalize the pivot column
//        normalize_pivot_col<<<nprGrid, linearBlock>>>(matrix_dev, result_dev, index, size);
		cudaDeviceSynchronize();

#ifdef DEBUG
//uint8_t matrix_host[size * size];
cudaMemcpy(matrix_host, matrix_dev, size * size, cudaMemcpyDeviceToHost);
printf("Step: normalize pivot row\n");
printf("matrix:\n");
show_squre_matrix_debug(matrix_host, size);
//uint8_t result_host[size * size];
cudaMemcpy(result_host, result_dev, size * size, cudaMemcpyDeviceToHost);
printf("result:\n");
show_squre_matrix_debug(result_host, size);
#endif

		dim3 ebrGrid(size, (int) (ceil((float) size / SINGLE_BLOCK_SIZE)));
		dim3 ebrBlock(1, min(size, SINGLE_BLOCK_SIZE)); 
        eliminate_by_row<<<ebrGrid, ebrBlock>>>(matrix_dev, result_dev, row, size);
		cudaDeviceSynchronize();

#ifdef DEBUG
//uint8_t matrix_host[size * size];
cudaMemcpy(matrix_host, matrix_dev, size * size, cudaMemcpyDeviceToHost);
printf("Step: eliminate by row\n");
printf("matrix:\n");
show_squre_matrix_debug(matrix_host, size);
//uint8_t result_host[size * size];
cudaMemcpy(result_host, result_dev, size * size, cudaMemcpyDeviceToHost);
printf("result:\n");
show_squre_matrix_debug(result_host, size);
#endif
    }

}

__global__ void gen_encoding_matrix(uint8_t *encodingMatrix, int row, int col)
{
	int i = threadIdx.x;
	int j = threadIdx.y;
	setup_tables();
	__syncthreads();
	encodingMatrix[i * col + j] = gf_pow((j + 1) % field_size, i);
}

__host__ float encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	int threadsPerBlock = 128;
	int tileWidthRow = parityBlockNum;
	int tileWidthCol = threadsPerBlock / tileWidthRow;
	int tileDepth = nativeBlockNum;
	int gridDimX = min((int) (ceil((float) chunkSize / tileWidthCol)), SINGLE_GRID_SIZE);
	int gridDimY = (int)( ceil((float) parityBlockNum / tileWidthRow));
	dim3 grid(gridDimX, gridDimY);
	dim3 block(tileWidthCol, tileWidthRow);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	int sMemMaxSize = deviceProp.sharedMemPerBlock;
	int sMemMinSize = (tileWidthRow + tileWidthCol) * tileDepth * sizeof(uint8_t);
	int tunedSMemSize = 2048;
	/* 
	size_t sMemSize = tunedSMemSize;
	if (sMemMinSize > tunedSMemSize)
	{
		sMemSize = sMemMinSize;
	}
	*/
	size_t sMemSize = sMemMinSize;
	float stepTime = 0;
	cudaEvent_t stepStart, stepStop;
	// create event
	cudaEventCreate(&stepStart);
	cudaEventCreate(&stepStop);
	// record event
	cudaEventRecord(stepStart);
	matrix_mul<<<grid, block, sMemSize>>>(parityCoeff, dataChunk, codeChunk, parityBlockNum, nativeBlockNum, chunkSize, tileWidthRow, tileWidthCol, tileDepth);
//	matrix_mul<<<grid, block>>>(parityCoeff, dataChunk, codeChunk, parityBlockNum, nativeBlockNum, chunkSize, tileWidthRow, tileWidthCol, tileDepth);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	return stepTime;
}

__host__ float decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	int threadsPerBlock = 128;
	int tileWidthRow = nativeBlockNum;
	if (tileWidthRow > 8)
	{
		tileWidthRow = 8;
	}
	int tileWidthCol = threadsPerBlock / tileWidthRow;
	int tileDepth = nativeBlockNum;
	int gridDimX = min((int) (ceil((float) chunkSize / tileWidthCol)), SINGLE_GRID_SIZE);
	int gridDimY = (int) (ceil((float) nativeBlockNum / tileWidthRow));
	dim3 grid(gridDimX, gridDimY);
	dim3 block(tileWidthCol, tileWidthRow);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	int sMemMaxSize = deviceProp.sharedMemPerBlock;
	int sMemMinSize = (tileWidthRow + tileWidthCol) * tileDepth * sizeof(uint8_t);
	int tunedSMemSize = 2048;
	/* 
	size_t sMemSize = tunedSMemSize;
	if (sMemMinSize > tunedSMemSize)
	{
		sMemSize = sMemMinSize;
	}
	*/
	size_t sMemSize = sMemMinSize;
	float stepTime = 0;
	cudaEvent_t stepStart, stepStop;
	// create event
	cudaEventCreate(&stepStart);
	cudaEventCreate(&stepStop);
	// record event
	cudaEventRecord(stepStart);
	matrix_mul<<<grid, block, sMemSize>>>(parityCoeff, codeChunk, dataChunk, nativeBlockNum, nativeBlockNum, chunkSize, tileWidthRow, tileWidthCol, tileDepth);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	return stepTime;
}

