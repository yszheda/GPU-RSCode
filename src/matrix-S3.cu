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

const int width = 8;
const int field_size = 1 << 8;
// const int field_size = 1 << width;

__shared__ uint8_t gflog[256];
__shared__ uint8_t gfexp[256];

__device__ int setup_tables()
{
//	const int width = 8;
//	const int field_size = 1 << width;
	const unsigned int prim_poly = 0435;
   	int log;
	int exp = 1;
	// use int as book-keeping index instead of unsigned int
	for (log = 0; log < field_size - 1; log++) 
	{
		if(exp > field_size) break;
		gflog[exp] = (uint8_t) log;
		gfexp[log] = (uint8_t) exp;
		exp = exp << 1;
		if (exp & field_size) 
		{
			exp = exp ^ prim_poly;
		}
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
	if (a == 0 || b == 0)
	{
		return 0;
	}
//	sum_log = (gflog[a] + gflog[b]) % (field_size-1);
	sum_log = gflog[a] + gflog[b];
	if (sum_log >= field_size - 1)
	{	
		sum_log -= field_size - 1;
	}
	return gfexp[sum_log];
}

__host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
{
	int sum_log;
	if (a == 0 || b == 0)
	{
		return 0;
	}
//	sum_log = (gflog[a] + gflog[b]) % (field_size-1);
	sum_log = gflog[a] + gflog[b];
	if (sum_log >= field_size - 1)
	{	
		sum_log -= field_size - 1;
	}
	return gfexp[sum_log];
}

__host__ __device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b)
{
	uint8_t sum_log = 0;
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
	uint8_t sum_log = 0;
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
//	diff_log = (gflog[a] - gflog[b]) % (field_size-1);
	diff_log = gflog[a] - gflog[b];
	if (diff_log < 0)
	{	
		diff_log += field_size - 1;
	}
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
//	diff_log = (gflog[a] - gflog[b]) % (field_size-1);
	diff_log = gflog[a] - gflog[b];
	if (diff_log < 0)
	{	
		diff_log += field_size - 1;
	}
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

	setup_tables();
	__syncthreads();

	bx = blockIdx.x;
	do {
		row = by*tileWidthRow + ty;
		col = bx*tileWidthCol + tx;
		product = 0;
		__syncthreads();
		
		if(row < n && col < m)
		{
//			for(int j = tx; j < tileDepth; j += blockDim.x)
//			{
//				sMem[ index(ty, j, tileDepth) ] = A[row*p + j];
//			}
			for(int j = ty; j < tileDepth; j += blockDim.y)
			{
				sMem[rowVectorSize + index(j, tx, tileWidthCol)] = B[col + j*m];
			}
			// Since blockDim.x > tileDepth for our applications,
			// we can fully parallelize loading matrix A into sMem.
			if (tx < tileDepth)
			{
				sMem[ index(ty, tx, tileDepth) ] = A[row*p + tx];
			}
//			if (ty < tileDepth)
//			{
//				sMem[rowVectorSize + index(ty, tx, tileWidthCol)] = B[col + ty*m];
//			}
			__syncthreads();
			
			for(int j = 0; j < tileDepth; j++)
			{
				product ^= gf_mul_bit(sMem[ index(ty, j, tileDepth) ], sMem[rowVectorSize + index(j, tx, tileWidthCol)]);
			}
			__syncthreads();
			C[row*m+col] = product;
		}
		bx += gridDim.x;
		col = bx*tileWidthCol + tx;
		__syncthreads();
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

	setup_tables();
	__syncthreads();

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

	setup_tables();
	__syncthreads();

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

	setup_tables();
	__syncthreads();

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

	setup_tables();
	__syncthreads();

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
	setup_tables();
	__syncthreads();
	encodingMatrix[i*col + j] = gf_pow((j+1) % field_size, i);
}

__host__ float encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	int threadsPerBlock = 128;
	int tileWidthRow = ceil((float)threadsPerBlock * parityBlockNum / (parityBlockNum + chunkSize));
	int tileWidthCol = threadsPerBlock / tileWidthRow;
	int tileDepth = nativeBlockNum;
	int gridDimX = min( (int)( ceil((float)chunkSize / tileWidthCol) ), SINGLE_GRID_SIZE );
	int gridDimY = (int)( ceil((float)parityBlockNum / tileWidthRow) );
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
	int tileWidthRow = ceil((float)threadsPerBlock * nativeBlockNum / (nativeBlockNum + chunkSize));
	int tileWidthCol = threadsPerBlock / tileWidthRow;
	int tileDepth = nativeBlockNum;
	int gridDimX = min( (int)( ceil((float)chunkSize / tileWidthCol) ), SINGLE_GRID_SIZE );
	int gridDimY = (int)( ceil((float)nativeBlockNum / tileWidthRow) );
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
//	matrix_mul<<<grid, block>>>(parityCoeff, codeChunk, dataChunk, nativeBlockNum, nativeBlockNum, chunkSize, tileWidthRow, tileWidthCol, tileDepth);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	return stepTime;
}

