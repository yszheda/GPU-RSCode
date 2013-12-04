#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <stdint.h>
//#include "galoisfield.h"

#define W 8
#define NW (1 << W) /* In other words, NW equals 2 to the w-th power */

//need tuning!!!
#define TILE_WIDTH 4
#define TILE_DEPTH 4

#define BUFFER_SIZE 256

#define DEBUG 

//#define BLOCKSIZE 16
//#define BLOCKSIZEMINUS1 15
#define BLOCKSIZE 4
#define BLOCKSIZEMINUS1 3

//#define USELOOPUNROLLING 1  
#define AVOIDBANKCONFLICTS 0    //this just runs faster :X

__shared__ uint8_t gflog[256];
__shared__ uint8_t gfexp[256];

//__global__ int setup_tables(int w, uint8_t *gflog, uint8_t *gfexp)
//__device__ int setup_tables(int w, uint8_t *gflog, uint8_t *gfexp)
__device__ int setup_tables(int w)
{
	unsigned int b;
   	unsigned int log;
	unsigned int x_to_w;
	unsigned int prim_poly;
//	unsigned int r;
//	unsigned int x;
//	unsigned int y;

	unsigned int prim_poly_4 = 023;
	unsigned int prim_poly_8 = 0435;
	//uint8_t prim_poly_8 = 285;
	unsigned int prim_poly_16 = 0210013;
	switch(w) 
	{
		case 4: prim_poly = prim_poly_4; break;
		case 8: prim_poly = prim_poly_8; break;
		case 16: prim_poly = prim_poly_16; break;
		default: return -1;
	}
	x_to_w = 1 << w;
	b = 1;
//	r = 0;
	for (log = 0; log < x_to_w-1; log++) 
	{
		/*
		r = 0;
		x = 1;
		y = log;
		while(y)
		{
			printf("y=%d\n",y);
			if(y & 1)
			{
				r = r ^ b;
			}
			y = y >> 1;
			x = x << 1;
			if (x & x_to_w) x = x ^ prim_poly;
		}
			printf("r=%d\n",r);
			printf("log=%d\n",log);
		*/
		if(b>x_to_w) break;
		gflog[b] = (uint8_t) log;
		gfexp[log] = (uint8_t) b;
		b = b << 1;
		if (b & x_to_w) 
		{
			b = b ^ prim_poly;
		}
	}
	return 0;
}

__device__ uint8_t gf_add(uint8_t a, uint8_t b)
{
	return a^b;
}

__device__ uint8_t gf_sub(uint8_t a, uint8_t b)
{
	return gf_add(a, b);
}

__device__ uint8_t gf_mul(uint8_t a, uint8_t b)
{
	uint8_t sum_log;
	if (a == 0 || b == 0)
	{
		return 0;
	}
	sum_log = (gflog[a] + gflog[b]) % (NW-1);
	/*
	sum_log = gflog[a] + gflog[b];
	if (sum_log >= NW-1)
	{	
		sum_log -= NW-1;
	}
	*/
	return gfexp[sum_log];
}

__device__ uint8_t gf_mul(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
{
	int sum_log;
	if (a == 0 || b == 0)
	{
		return 0;
	}
	sum_log = (gflog[a] + gflog[b]) % (NW-1);
	/*
	sum_log = gflog[a] + gflog[b];
	if (sum_log >= NW-1)
	{	
		sum_log -= NW-1;
	}
	*/
	return gfexp[sum_log];
}

uint8_t gf_mul_bit(uint8_t a, uint8_t b)
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

__device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
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

__device__ uint8_t gf_div(uint8_t a, uint8_t b)
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
//	diff_log = (gflog[a] - gflog[b]) % (NW-1);
	diff_log = gflog[a] - gflog[b];
	if (diff_log < 0)
	{	
		diff_log += NW-1;
	}
	return gfexp[diff_log];
}

__device__ uint8_t gf_div(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
{
	uint8_t diff_log;
	if (a == 0)
	{	
		return 0;
	}
	/* Can’t divide by 0 */
	if (b == 0)
	{
		return -1;
	}
	diff_log = (gflog[a] - gflog[b]) % (NW-1);
	/*
	diff_log = gflog[a] - gflog[b];
	if (diff_log < 0)
	{	
		diff_log += NW-1;
	}
	*/
	return gfexp[diff_log];
}

__device__ uint8_t gf_pow(uint8_t a, uint8_t power)
{
	uint8_t pow_log = (gflog[a] * power) % (NW-1);
	return gfexp[pow_log];
}

__device__ uint8_t gf_pow(uint8_t a, uint8_t power, uint8_t *gflog, uint8_t *gfexp)
{
	uint8_t pow_log = (gflog[a] * power) % (NW-1);
	return gfexp[pow_log];
}


// C=AB
// A: nxp
// B: pxm
// C: nxm
__device__ void matrix_mul(unsigned char *A, unsigned char *B, unsigned char *C, int n, int p, int m)
{
	__shared__ int rowVector[TILE_WIDTH][TILE_DEPTH];
	__shared__ int colVector[TILE_DEPTH][TILE_WIDTH];
	__shared__ int product[TILE_WIDTH][TILE_WIDTH];

	int bx = blockIdx.x;
   	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int row;
	int col;
	int px;
	int py;	

	setup_tables(8);
	__syncthreads();

	for(py=ty; py<TILE_WIDTH; py+=blockDim.y)
	{
		for(px=tx; px<TILE_WIDTH; px+=blockDim.x)
		{
			row = by*TILE_WIDTH+py;
			col = bx*TILE_WIDTH+px;
			product[py][px] = 0;
			__syncthreads();
		
			for(int i=0; i<(int)(ceil((float)p/TILE_DEPTH)); i++)
			{
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
				__syncthreads();
			}
			C[row*m+col] = product[py][px];
		}
	}
	/*
	int i;
	int j;
	int k;
	setup_tables(8);
	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
		{
			for(k=0; k<p; k++)
			{
				C[i*m+j] = gf_add( C[i*m+j], gf_mul( A[i*p+k], B[k*m+j] ) );
			}
		}
	}
	*/
}

__global__ void decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	matrix_mul(parityCoeff, codeChunk, dataChunk, nativeBlockNum, nativeBlockNum, chunkSize);
}

/********************************************************************
*  File: GPUadjustRow_kernel.cu
*
*  Description: 
*   adjust the rest of row with the pre-calculated pivot elements from step 1.
*   Also adjust the inverse block with the identity matrix.
*	
*  Includes 2 Functions:
*    - adjustColL_kernel
*    - adjustColU_kernel   
*  
*
*  Arguments: 
*	- float *dInData      Input Matrix 1D, on the device
*   - int size            Matrix dimension in size
*	
* Used custom Routines:
*	- 
*
*********************************************************************/
//#include <stdio.h>
//#include <cuda.h>
//#include <stdint.h>
//
//#ifndef GPUGAUSSEIDEL_H
//#define GPUGAUSSEIDEL_H
//
//#ifdef __cplusplus
//   extern "C" {
//#endif
//
//#define BLOCKSIZE 16
//#define BLOCKSIZEMINUS1 15
//
//#define USELOOPUNROLLING 1  
//#define AVOIDBANKCONFLICTS 0    //this just runs faster :X
//
//int GPUGausSeidel (float* matrix, 
//                   float* output, 
//                   int size);

//#ifdef __cplusplus
//   }
//#endif
//
//#endif


//#ifndef GPUADJUSTROW_KERNEL_H
//#define GPUADJUSTROW_KERNEL_H

__global__ void adjustRowL_kernel (uint8_t *dMatrixIn, uint8_t *dMatrixInDiag,
                                   uint8_t *dMatrixInv, int width, int diagEl)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int bx = blockIdx.x;

    __shared__ uint8_t pivotBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];
    __shared__ uint8_t inBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];
    __shared__ uint8_t invBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];


	setup_tables(8);
	__syncthreads();

    /*
     * Adjust the rest blocks which are right from the prepared block of step 1
     * and adjust the inverse blocks
     */
    if (bx * BLOCKSIZE > (diagEl + 1))
    {
        pivotBlock[ty][tx] = dMatrixInDiag[ty * width + tx];
        inBlock[ty][tx] = dMatrixIn[ty * width + bx * BLOCKSIZE +tx];
        invBlock[ty][tx] = dMatrixInv[ty * width + bx * BLOCKSIZE +tx];

        __syncthreads ();
#ifdef USELOOPUNROLLING
#pragma unroll BLOCKSIZE
#endif
        //i equals the current row where the pivot elements are stored
        for (int i = 0; i < BLOCKSIZEMINUS1; i++)
        {
            // if the cols are below  
            if (ty > i)
            {
                uint8_t pivot = pivotBlock[ty][i];
                //Subtract the row
                inBlock[ty][tx] ^= gf_mul(inBlock[i][tx], pivot);
                invBlock[ty][tx] ^= gf_mul(invBlock[i][tx], pivot);
            }

            __syncthreads ();
        }
        //Store the results back in device memory
        dMatrixIn[ty * width + bx * BLOCKSIZE +tx] = inBlock[ty][tx];
        dMatrixInv[ty * width + bx * BLOCKSIZE +tx] = invBlock[ty][tx];
    }
    /*
     * Adjust the last blocks from the indentity matrix which are left 
     */
    else
    {
        pivotBlock[ty][tx] = dMatrixInDiag[ty * width + tx];
        invBlock[ty][tx] = dMatrixInv[ty * width + bx * BLOCKSIZE +tx];

        __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif
        for (int i = 0; i < BLOCKSIZEMINUS1; i++)//last changed
        {
            if (ty > i)
            {
                uint8_t pivot = pivotBlock[ty][i];

                invBlock[ty][tx] ^= gf_mul(invBlock[i][tx], pivot);
            }

            __syncthreads ();
        }
        dMatrixInv[ty * width + bx * BLOCKSIZE + tx] = invBlock[ty][tx];
    }
}



__global__ void adjustRowU_kernel (uint8_t *dMatrixIn, uint8_t *dMatrixInv, int width,
                                   int diagEl)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int bx = blockIdx.x;

    __shared__ uint8_t pivotBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];
    __shared__ uint8_t invBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];


	setup_tables(8);
	__syncthreads();

    pivotBlock[ty][tx] = dMatrixIn[ty * width + tx];
    invBlock[ty][tx] = dMatrixInv[ty * width + bx * BLOCKSIZE +tx];

    __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif
    for (int i = BLOCKSIZEMINUS1; i > 0; i--)
    {
        if (ty < i)
        {
            uint8_t pivot = pivotBlock[ty][i];

            invBlock[ty][tx] ^= gf_mul(invBlock[i][tx], pivot);
        }

        __syncthreads ();
    }

    dMatrixInv[ty * width + bx * BLOCKSIZE +tx] = invBlock[ty][tx];
}

//#endif
/********************************************************************
*  File: GPUeliminateBlock_kernel.cu
*
*  Description: 
*   Calculate the Gaus-Seidel algorithm on a small BLOCK of the Inputmatrix. 
*   To benefit from the Shared Memory the BLOCK is sized to BLOCKSIZE = 16.
*	
*  Includes 2 Functions:
*    - eliminateBlockL_kernel
*    - eliminateBlockU_kernel   
*  
*
*  Arguments: 
*	- float *dInData      Input Matrix 1D, on the device
*   - int size            Matrix dimension in size
*	
* Used custom Routines:
*	- 
*
*********************************************************************/


//#ifndef ELIMINATEBLOCK_H
//#define ELIMINATEBLOCK_H


__global__ void eliminateBlockL_kernel (uint8_t *dInData, 
                                        int size)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    __shared__ uint8_t triangleBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];


	setup_tables(8);
	__syncthreads();

    triangleBlock[ty][tx] = dInData[ty * size + tx];
    __syncthreads ();

    
#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZEMINUS1
#endif
    //i equals the current row
    for (int i = 0; i < BLOCKSIZEMINUS1; i++)
    {
        // calculate the pivot element to get the current row i to zero
        uint8_t pivotEl = gf_div(triangleBlock[ty][i], triangleBlock[i][i]);

        __syncthreads ();       // Each pivotEl have to be calculated and store in the registers

        if (ty > i)             // If all cols (ty) are below the current row (i)?
        {
            if (tx > i)         // The element is right to the current row, subtract the element
            {
                triangleBlock[ty][tx] ^= gf_mul(pivotEl, triangleBlock[i][tx]);
            }
            if (tx == i)        // Store the pivot element in the current row
            {
                triangleBlock[ty][tx] = pivotEl;
            }
        }
        __syncthreads ();       // Wait for each thread
    }

    dInData[ty * size + tx] = triangleBlock[ty][tx];       // Write the result back to memory
}


__global__ void eliminateBlockU_kernel (uint8_t *dInData, int size)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    __shared__ uint8_t triangleBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];


	setup_tables(8);
	__syncthreads();

    triangleBlock[ty][tx] = dInData[ty * size + tx];
    __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZEMINUS1
#endif    
    //i equals the current row
    for (int i = BLOCKSIZEMINUS1; i > 0; i--)
    {
        // calculate the pivot element to get the current row i to zero
        uint8_t pivotEl = gf_div(triangleBlock[ty][i], triangleBlock[i][i]);

        __syncthreads ();       // Each pivotEl have to be calculated and store in the registers

        if (ty < i)             // If all rows (ty) are above the current row (i)?
        {
            if (tx < i)         // The element is left to the current row, subtract the element
            {
                triangleBlock[ty][tx] ^= gf_mul(pivotEl, triangleBlock[i][tx]);
            }
            if (tx == i)        // Store the pivot element in the current row
            {
                triangleBlock[ty][tx] = pivotEl;
            }
        }
        __syncthreads ();        // Wait for each thread
    }

    dInData[ty * size + tx] = triangleBlock[ty][tx];       //Write the result back to device memory
}

//#endif
/********************************************************************
*  File: GPUeliminateCol_kernel.cu
*
*  Description: Calculate each pivot elements above/over (L/U) the block
*               from step 1. And fill The blocks with them
*   
*	
*  Includes 2 Functions:
*    - eliminateColL_kernel
*    - eliminateColU_kernel   
*  
*
*  Arguments: 
*	- float *dMatrixIn      Input Matrix 1D pointed to the current row, on the device 
*   - int size              Matrix dimension in size
*   - int diagEl            The adjusted blockoffset from step 1
*	
* Used custom Routines:
*	- 
*
*********************************************************************/
//#ifndef ELIMINATE_COL_KERNEL
//#define ELIMINATE_COL_KERNEL


__global__ void eliminateColL_kernel (uint8_t *dMatrixIn, int size, int diagEl)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    //bx is used to adress the Blocks above the precalculated block from step 1
    int bx = blockIdx.x;    


	setup_tables(8);
	__syncthreads();

    //only the blocks can enter which are above the precalculated block from step 1
    if (bx * BLOCKSIZE > (diagEl + 1))
    {
        int offset = diagEl * size;
        int blockOffset = bx * BLOCKSIZE *size;

        __shared__ uint8_t pivotBlock[BLOCKSIZE][BLOCKSIZE];
        __shared__ uint8_t inBlock[BLOCKSIZE][BLOCKSIZE];

        pivotBlock[ty][tx] = dMatrixIn[offset + ty * size + tx];   // The Block from step 1
        inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx]; // each Block which is above the pivotBlock

        __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif
        //iterate through the block und calculate the pivot elements
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            uint8_t pivotEl = gf_div(inBlock[ty][i], pivotBlock[i][i]);

            __syncthreads ();

            //adjust all values right to the current interation step
            if (tx > i)
            {
                //substract the row
                inBlock[ty][tx] ^= gf_mul(pivotBlock[i][tx], pivotEl);
            }
            //store the pivot element in the col
            else
            {
                inBlock[ty][i] = pivotEl;
            }

            __syncthreads ();
        }

        dMatrixIn[blockOffset + ty * size + tx] = inBlock[ty][tx];
    }
}


__global__ void eliminateColU_kernel (uint8_t *dMatrixIn, int size, int diagEl)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int bx = blockIdx.x;


	setup_tables(8);
	__syncthreads();

    if (bx * BLOCKSIZE < diagEl)
    {
        int offset = diagEl * size;
        int blockOffset = bx * BLOCKSIZE *size;

        __shared__ uint8_t pivotBlock[BLOCKSIZE][BLOCKSIZE];
        __shared__ uint8_t inBlock[BLOCKSIZE][BLOCKSIZE];

        pivotBlock[ty][tx] = dMatrixIn[offset + ty * size + tx];
        inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx];

        __syncthreads ();

#ifdef USELOOPUNROLLING
#pragma unroll BLOCKSIZE
#endif

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            uint8_t pivotEl = gf_div(inBlock[ty][i], pivotBlock[i][i]);

            __syncthreads ();

            if (tx < i)
            {
                inBlock[ty][tx] ^= gf_mul(pivotBlock[i][tx], pivotEl);
            }
            else/* if (tx == i)*/
            {
                inBlock[ty][i] = pivotEl;
            }

            __syncthreads ();
        }

        dMatrixIn[blockOffset + ty * size + tx] = inBlock[ty][tx];
    }
}

//#endif
/********************************************************************
*  File: GPUeliminateRest_kernel.cu
*
*  Description: Adjust the rest of the entire inv-/matrix with the precalculated
*               pivot elements from step 3. 
*   
*	
*  Includes 2 Functions:
*    - eliminateRestL_kernel
*    - eliminateRestU_kernel   
*  
*
*  Arguments: 
*	- float *dMatrixIn      Input Matrix 1D 
*	- float *dMatrixInv     Invers Matrix 1D
*   - int size              Matrix dimension in size
*   - int diagEl            The adjusted blockoffset from step 1
*	
* Used custom Routines:
*	- 
*
*********************************************************************/
//#ifndef ELIMINATE_REST_KERNEL
//#define ELIMINATE_REST_KERNEL


__global__ void eliminateRestL_kernel (uint8_t *dMatrixIn, uint8_t *dMatrixInv, int size,
                                       int diagEl)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int bx = blockIdx.x;
    int by = blockIdx.y;

    __shared__ uint8_t pivEl[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];
    __shared__ uint8_t pivBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];
    __shared__ uint8_t inBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];


	setup_tables(8);
	__syncthreads();

    //rest of the unadjusted Matrix which is right above the diagEl
    if (bx * BLOCKSIZE > (diagEl + 1) && by * BLOCKSIZE > (diagEl + 1))
    {
        int blockOffset = by * BLOCKSIZE * size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE * size + diagEl;
        int blockPivOffset = diagEl * size + bx * BLOCKSIZE;

        inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx];
        pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixIn[blockPivOffset + ty * size + tx];
        __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif
        //Subtract each row from the input Matrix =>dMatrixIn
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            inBlock[ty][tx] ^= gf_mul(pivEl[ty][i], pivBlock[i][tx]);
        }
        __syncthreads ();
        dMatrixIn[blockOffset + ty * size + tx] = inBlock[ty][tx];
        __syncthreads ();

        inBlock[ty][tx] = dMatrixInv[blockOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixInv[blockPivOffset + ty * size + tx];

        __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif
        //Subtract each row from the invers Matrix =>dMatrixInv
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            inBlock[ty][tx] ^= gf_mul(pivEl[ty][i], pivBlock[i][tx]);
        }

        __syncthreads ();
        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
    }
    //Adjust the left Blocks from the invers matrix which are left from the diagEl
    else if (by * BLOCKSIZE > (diagEl + 1))
    {
        int blockOffset = by * BLOCKSIZE * size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE * size + diagEl;
        int blockPivOffset = diagEl * size + bx * BLOCKSIZE;

        pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        inBlock[ty][tx] = dMatrixInv[blockOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixInv[blockPivOffset + ty * size + tx];
        __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif

        for (int i = 0; i < BLOCKSIZE; i++)
        {
            inBlock[ty][tx] ^= gf_mul(pivEl[ty][i], pivBlock[i][tx]);
        }

        __syncthreads ();
        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
    }
}


__global__ void eliminateRestU_kernel (uint8_t *dMatrixIn, uint8_t *dMatrixInv, int size,
                                       int diagEl)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int bx = blockIdx.x;
    int by = blockIdx.y;

    __shared__ uint8_t pivEl[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];
    __shared__ uint8_t pivBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];
    __shared__ uint8_t inBlock[BLOCKSIZE][BLOCKSIZE+AVOIDBANKCONFLICTS];


	setup_tables(8);
	__syncthreads();

    //rest der unbearbeiteten Matrix bearbeiten
    if ((bx * BLOCKSIZE + 1) <diagEl && (by * BLOCKSIZE +1) <diagEl)     //linke seite von in; 0-pivblock
    {

        int blockOffset = by * BLOCKSIZE * size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE * size + diagEl;
        int blockPivOffset = diagEl * size + bx * BLOCKSIZE;

        inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx];
        pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixIn[blockPivOffset + ty * size + tx];
        __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            inBlock[ty][tx] ^= gf_mul(pivEl[ty][i], pivBlock[i][tx]);
        }
        __syncthreads ();
        dMatrixIn[blockOffset + ty * size + tx] = inBlock[ty][tx];
        __syncthreads ();

        inBlock[ty][tx] = dMatrixInv[blockOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixInv[blockPivOffset + ty * size + tx];

        __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            inBlock[ty][tx] ^= gf_mul(pivEl[ty][i], pivBlock[i][tx]);
        }

        __syncthreads ();
        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
    }
    else if (by * BLOCKSIZE <(diagEl))
    {
        int blockOffset = by * BLOCKSIZE *size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE *size + diagEl;
        int blockPivOffset = diagEl * size + bx * BLOCKSIZE;


        pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        inBlock[ty][tx] = dMatrixInv[blockOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixInv[blockPivOffset + ty * size + tx];
        __syncthreads ();

#ifdef USELOOPUNROLLING
    #pragma unroll BLOCKSIZE
#endif

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            inBlock[ty][tx] ^= gf_mul(pivEl[ty][i], pivBlock[i][tx]);
        }

        __syncthreads ();
        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
    }
}

//#endif
/********************************************************************
*  File: GPUnormalizeDiag.cu
*
*  Description: Force the entries of diagonal matrix to 1 to get a
*               left sided entity matrix.
*   
*	
*  Includes 1 Functions:
*    - normalizeDiag_kernel 
*  
*
*  Arguments: 
*	- float *diagMatrix    Input Matrix 1D 
*	- float *invMatrix     Invers Matrix 1D
*   - int size             Matrix dimension in size
*   - int row              The current Block offset of the row
*	
* Used custom Routines:
*	- 
*
*********************************************************************/
//#ifndef NORMALIZEDIAG_CU
//#define NORMALIZEDIAG_CU

__global__ void normalizeDiag_kernel (uint8_t *diagMatrix, uint8_t *invMatrix, int size,
                                      int row)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int bx = blockIdx.x;

    int blockOffset = bx * BLOCKSIZE;
    __shared__ uint8_t diagEl[BLOCKSIZE];


	setup_tables(8);
	__syncthreads();

    if (tx == ty)
    {
        diagEl[ty] = diagMatrix[row + ty * size + tx];
    }
    __syncthreads ();

    invMatrix[blockOffset + ty * size + tx] =
                            gf_div(invMatrix[blockOffset + ty * size + tx], diagEl[ty]);
}

//#endif
/********************************************************************
*  File: GPUnormalizeDiag.cu
*
*  Description: Generate an identitymatrix
*   
*	
*  Includes 1 Functions:
*    - GPUsetIdentity 
*  
*
*  Arguments: 
*	- float *diagMatrix    Input Matrix 1D 
*	- float *invMatrix     Invers Matrix 1D
*   - int size             Matrix dimension in size
*   - int row              The current Block offset of the row
*	
* Used custom Routines:
*	- 
*
*********************************************************************/
//#ifndef GPUSETIDENTITY_KERNEL_H
//#define GPUSETIDENTITY_KERNEL_H

__global__ void GPUsetIdentity (uint8_t* matrix,
                                int width)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;

    int offset = bx * BLOCKSIZE + tx;
    matrix[offset * width + offset] = 1;
}

//#endif
/********************************************************************
*  File: GPUGausSeidel.c
*
*  Description:
*	
*	Mainfunction to compute an inverse Matrix from a positive definite 
*   Matrix on the GPU. The Routine is using the Gaus Seidel Matrix invertion 
*   algorithm. 
*   
*   1   2   1   |  1  0   0               1   0   0  |  -2.5   1.5   0.5
*   2   3   1   |  0  1   0       =>      0   1   0  |   1.5  -0.5  -0.5 
*   1   1   2   |  0  0   1               0   0   1  |   0.5  -0.5   0.5 
*   Inputmatrix       E                       E          Inverse Matrix
*
*  Arguments: 
*	- float *hDataIn      Input Matrix 1D, no data changes
*   - float *hDataOut     Output Matrix 1D, the inverse datamatrix  
*   - int size            Matrix dimension in size, width = height = size
*	
*  Used custom kernels rutines:
*	- GPUsetIdentity          
*   - eliminateBlockL_kernel    
*   - adjustColL_kernel         
*   - eliminateColL_kernel      
*   - eliminateRestL_kernel     
*
*   - eliminateBlockU_kernel    
*   - adjustColU_kernel         
*   - eliminateColU_kernel      
*   - eliminateRestU_kernel     
*
*********************************************************************/

//#ifndef GPUGAUSSEIDEL_H
//#define GPUGAUSSEIDEL_H
//
//#ifdef __cplusplus
//   extern "C" {
//#endif
//
//#define BLOCKSIZE 16
//#define BLOCKSIZEMINUS1 15
//
//#define USELOOPUNROLLING 1  
//#define AVOIDBANKCONFLICTS 0    //this just runs faster :X
//
//int GPUGausSeidel (float* matrix, 
//                   float* output, 
//                   int size);
//
//#ifdef __cplusplus
//   }
//#endif
//
//#endif

//#include <cutil.h>
//
//#include "../Header/GPUGausSeidel.h"
//
//#include "GPUsetIdentity_kernel.cu"
//#include "GPUeliminateBlock_kernel.cu"
//#include "GPUeliminateCol_kernel.cu"
//#include "GPUeliminateRest_kernel.cu"
//#include "GPUadjustRow_kernel.cu"
//#include "GPUnormalizeDiag_kernel.cu"


void GPUGausSeidel (uint8_t *dDataIn, uint8_t *dDataInv, int size)
//void GPUGausSeidel (uint8_t *hDataIn, uint8_t *hDataOut)
{
    int i;
/*
    uint8_t *dDataIn;
    uint8_t *dDataInv;
    int size2InBytes = size * size * sizeof (uint8_t);
  
    //Allocating memory for the datamatrix and identity matrix (Einheitsmatrix)
    cudaMalloc ((void **) &dDataIn, size2InBytes);
    cudaMalloc ((void **) &dDataInv, size2InBytes);
    
    //Prepare the calculation of the identitymatrix
    cudaMemset ((void *) dDataInv, 0, size2InBytes);
    //Transfair the matrix from host to device
    cudaMemcpy ((void *) dDataIn, (void *) hDataIn, size2InBytes,
                cudaMemcpyHostToDevice);
*/
    //Used SP/MP for calculations
    dim3 idyThreads (BLOCKSIZE);    
    dim3 idyBlocks (size / BLOCKSIZE);
    dim3 nThreads (BLOCKSIZE, BLOCKSIZE);   
    dim3 nBlocks (size / BLOCKSIZE);        
    dim3 nBlocksRest (size / BLOCKSIZE, size / BLOCKSIZE);

    //Calculate the Identitymatrix 
    GPUsetIdentity <<< idyBlocks, idyThreads >>> (dDataInv, size);
    cudaThreadSynchronize ();
	cudaDeviceSynchronize();

    //calculate the right diagonal Matrix (L)
    for (i = 0; i < size; i += BLOCKSIZE)
    {
        int offset = i * size + i;

        /* step 1:
         *  calculate the triangle matrix
         *  store the pivot elements to left part of the triangel
         */

        eliminateBlockL_kernel <<< 1, nThreads >>> (dDataIn + offset, size);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();

        /* step 2:
         *  calculate the rest of the rows with the pivot elements from step 1
         *  
         */
        adjustRowL_kernel <<< nBlocks, nThreads >>> (dDataIn + i * size, dDataIn + offset,
                                                     dDataInv + i * size, size, i);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();


        /* step 3:
         *Fill the colls below the block with the pivot elements they are used
         *    to get the colls to zero and multiply with the row
         */
        eliminateColL_kernel <<< nBlocks, nThreads >>> (dDataIn + i, size, i);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();

        /* step 4:
         *  Adjust the rest of the Matrix with the calculated pivot Elements
         *  El_new_0 -= (p0+p1+p2..+p15) * El_piv_0
         */
        eliminateRestL_kernel <<< nBlocksRest, nThreads >>> (dDataIn, dDataInv, size, i);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();
    }

    //Set the left lower diagonalmatrix to zero (async?)
    for (i = 1; i < size; i++)
    {
        int offset = i * size;
        cudaMemset ((void *) (dDataIn + offset), 0, i * sizeof (uint8_t));
    }
    cudaThreadSynchronize ();
	cudaDeviceSynchronize();


    //calculate the right diagonal Matrix (U)
    for (i = (size - BLOCKSIZE); i >= 0; i -= BLOCKSIZE)
    {
        int offset = i * size + i;

        /* step 1:
         *  calculate the triangle matrix
         *  store the pivot elements to left part of the triangel
         */
        eliminateBlockU_kernel <<< 1, nThreads >>> (dDataIn + offset, size);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();

        /* step 2:
         *  calculate the rest of the rows with the pivot elements from step 1
         *  
         */
        adjustRowU_kernel <<< nBlocks, nThreads >>> (dDataIn + offset,
                                                     dDataInv + i * size, size, i);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();

        /* step 3:
         *  Fill the colls below the block with the pivot elements they are used
         *      to get the colls to zero and multiply with the row
         */
        eliminateColU_kernel <<< nBlocks, nThreads >>> (dDataIn + i, size, i);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();

        /* step 4:
         *  Adjust the rest of the Matrix with the calculated pivot Elements
         *  El_new_0 -= (p0+p1+p2..+p15) * El_piv_0
         */
        eliminateRestU_kernel <<< nBlocksRest, nThreads >>> (dDataIn, dDataInv, size, i);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();
    }
    
    /*
     * force the diagonal entries to 1
     */
    for (i = 0; i < size; i += BLOCKSIZE)
    {
        int rowOffset = i * size;
        normalizeDiag_kernel <<< nBlocks, nThreads >>> (dDataIn + rowOffset,
                                                        dDataInv + rowOffset, size, i);
        cudaThreadSynchronize ();
	cudaDeviceSynchronize();
    }

/*
    cudaMemcpy ((void *) hDataOut, (void *) dDataInv, size2InBytes,
                cudaMemcpyDeviceToHost);

    cudaFree (dDataIn);
    cudaFree (dDataInv);

    return 0;
*/
}

void show_decoding_matrix(uint8_t *decodingMatrix, int size)
{
	int i;
	int j;
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			printf("%d ", decodingMatrix[i*size+j]);
		}
		printf("\n");
	}
}
int main()
{
	int nativeBlockNum = 4;
	int parityBlockNum = 2;
	int chunkSize = 1;

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host
	uint8_t *dataBuf_d;		//device
	uint8_t *codeBuf_d;		//device

	int dataSize;
	int codeSize;

	FILE *fp_in;
	FILE *fp_out;
	if( ( fp_in = fopen("native_0","rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}
	fseek(fp_in, 0L, SEEK_END);
	chunkSize = ftell(fp_in);


	dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( dataSize );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( codeSize);
	memset(codeBuf, 0, codeSize);
	cudaMalloc( (void **)&dataBuf_d, dataSize );
	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc( (void **)&codeBuf_d, codeSize );
	cudaMemset(codeBuf_d, 0, codeSize);

/*	
	int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	int codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( dataSize );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( codeSize);
	memset(codeBuf, 0, codeSize);
	cudaMalloc( (void **)&dataBuf_d, dataSize );
	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc( (void **)&codeBuf_d, codeSize );
	cudaMemset(codeBuf_d, 0, codeSize);

//
	FILE *fp_in;
	FILE *fp_out;
	if( ( fp_in = fopen(argv[1],"rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}
	fseek(fp_in, 0L, SEEK_END);
	//ftell() get the total size of the file
	chunkSize = (int) (ceil( (float)ftell(fp_in) / nativeBlockNum )); 
//	
	cudaMemcpy(codeBuf_d, codeBuf, codeSize, cudaMemcpyHostToDevice);
//
*/
#ifdef DEBUG
//	int matrixSize = nativeBlockNum*nativeBlockNum*sizeof(uint8_t);
	int matrixSize = 16;
	uint8_t testMatrix[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {1, 1, 1, 1}, {1, 2, 3, 4}};
//	uint8_t testMatrix[4][4] = {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 1, 1, 1}, {1, 2, 3, 4}};
//	uint8_t testMatrix[16] = {0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 3, 4};
//	uint8_t testMatrix[16] = {1, 1, 1, 1, 1, 2, 3, 4, 0, 0, 1, 0, 0, 0, 0, 1};
	uint8_t *encodingMatrix_d;	//device
	uint8_t *decodingMatrix;	//host
	uint8_t *decodingMatrix_d;	//device
	decodingMatrix = (uint8_t*) malloc( matrixSize );
	cudaMalloc( (void **)&encodingMatrix_d, matrixSize );
	cudaMalloc( (void **)&decodingMatrix_d, matrixSize );
	cudaMemcpy(encodingMatrix_d, testMatrix, matrixSize, cudaMemcpyHostToDevice);
	GPUGausSeidel (encodingMatrix_d, decodingMatrix_d, nativeBlockNum);


//    int i;
//	int size = nativeBlockNum;
///*
//    uint8_t *dDataIn;
//    uint8_t *dDataInv;
//    int size2InBytes = size * size * sizeof (uint8_t);
//  
//    //Allocating memory for the datamatrix and identity matrix (Einheitsmatrix)
//    cudaMalloc ((void **) &dDataIn, size2InBytes);
//    cudaMalloc ((void **) &dDataInv, size2InBytes);
//    
//    //Prepare the calculation of the identitymatrix
//    cudaMemset ((void *) dDataInv, 0, size2InBytes);
//    //Transfair the matrix from host to device
//    cudaMemcpy ((void *) dDataIn, (void *) hDataIn, size2InBytes,
//                cudaMemcpyHostToDevice);
//*/
//    //Used SP/MP for calculations
//    dim3 idyThreads (BLOCKSIZE);    
//    dim3 idyBlocks (size / BLOCKSIZE);
//    dim3 nThreads (BLOCKSIZE, BLOCKSIZE);   
//    dim3 nBlocks (size / BLOCKSIZE);        
//    dim3 nBlocksRest (size / BLOCKSIZE, size / BLOCKSIZE);
//
//    //Calculate the Identitymatrix 
//    GPUsetIdentity <<< idyBlocks, idyThreads >>> (decodingMatrix_d, size);
//    cudaThreadSynchronize ();
//
//    //calculate the right diagonal Matrix (L)
//    for (i = 0; i < size; i += BLOCKSIZE)
//    {
//        int offset = i * size + i;
//
//        /* step 1:
//         *  calculate the triangle matrix
//         *  store the pivot elements to left part of the triangel
//         */
//
//        eliminateBlockL_kernel <<< 1, nThreads >>> (encodingMatrix_d + offset, size);
//        cudaThreadSynchronize ();
//
//        /* step 2:
//         *  calculate the rest of the rows with the pivot elements from step 1
//         *  
//         */
//        adjustRowL_kernel <<< nBlocks, nThreads >>> (encodingMatrix_d + i * size, encodingMatrix_d + offset,
//                                                     decodingMatrix_d + i * size, size, i);
//        cudaThreadSynchronize ();
//
//
//        /* step 3:
//         *Fill the colls below the block with the pivot elements they are used
//         *    to get the colls to zero and multiply with the row
//         */
//        eliminateColL_kernel <<< nBlocks, nThreads >>> (decodingMatrix_d + i, size, i);
//        cudaThreadSynchronize ();
//
//        /* step 4:
//         *  Adjust the rest of the Matrix with the calculated pivot Elements
//         *  El_new_0 -= (p0+p1+p2..+p15) * El_piv_0
//         */
//        eliminateRestL_kernel <<< nBlocksRest, nThreads >>> (encodingMatrix_d, decodingMatrix_d, size, i);
//        cudaThreadSynchronize ();
//    }
//
//    //Set the left lower diagonalmatrix to zero (async?)
//    for (i = 1; i < size; i++)
//    {
//        int offset = i * size;
//        cudaMemset ((void *) (encodingMatrix_d + offset), 0, i * sizeof (uint8_t));
//    }
//    cudaThreadSynchronize ();
//
//
//    //calculate the right diagonal Matrix (U)
//    for (i = (size - BLOCKSIZE); i >= 0; i -= BLOCKSIZE)
//    {
//        int offset = i * size + i;
//
//        /* step 1:
//         *  calculate the triangle matrix
//         *  store the pivot elements to left part of the triangel
//         */
//        eliminateBlockU_kernel <<< 1, nThreads >>> (encodingMatrix_d + offset, size);
//        cudaThreadSynchronize ();
//
//        /* step 2:
//         *  calculate the rest of the rows with the pivot elements from step 1
//         *  
//         */
//        adjustRowU_kernel <<< nBlocks, nThreads >>> (encodingMatrix_d + offset,
//                                                     decodingMatrix_d + i * size, size, i);
//        cudaThreadSynchronize ();
//
//        /* step 3:
//         *  Fill the colls below the block with the pivot elements they are used
//         *      to get the colls to zero and multiply with the row
//         */
//        eliminateColU_kernel <<< nBlocks, nThreads >>> (encodingMatrix_d + i, size, i);
//        cudaThreadSynchronize ();
//
//        /* step 4:
//         *  Adjust the rest of the Matrix with the calculated pivot Elements
//         *  El_new_0 -= (p0+p1+p2..+p15) * El_piv_0
//         */
//        eliminateRestU_kernel <<< nBlocksRest, nThreads >>> (encodingMatrix_d, decodingMatrix_d, size, i);
//        cudaThreadSynchronize ();
//    }
//    
//    /*
//     * force the diagonal entries to 1
//     */
//    for (i = 0; i < size; i += BLOCKSIZE)
//    {
//        int rowOffset = i * size;
//        normalizeDiag_kernel <<< nBlocks, nThreads >>> (encodingMatrix_d + rowOffset,
//                                                        decodingMatrix_d + rowOffset, size, i);
//        cudaThreadSynchronize ();
//    }
//










	cudaMemcpy(decodingMatrix, decodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
	show_decoding_matrix(decodingMatrix, nativeBlockNum);
/*
	FILE *fp_in;
	FILE *fp_out;
	if( ( fp_in = fopen("native_2","rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}
	fseek(fp_in, 0L, SEEK_END);
	chunkSize = ftell(fp_in);


	dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( dataSize );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( codeSize);
	memset(codeBuf, 0, codeSize);
	cudaMalloc( (void **)&dataBuf_d, dataSize );
	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc( (void **)&codeBuf_d, codeSize );
	cudaMemset(codeBuf_d, 0, codeSize);
*/
//given a test matrix
//	uint8_t decodingMatrix_h[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {78, 104, 167, 180}, {79, 103, 166, 105}};
	uint8_t decodingMatrix_h[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {104, 187, 210, 186}, {105, 186, 211, 186}};
	cudaMemcpy(decodingMatrix_d, decodingMatrix_h, matrixSize, cudaMemcpyHostToDevice);
	cudaMemcpy(decodingMatrix, decodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
	show_decoding_matrix(decodingMatrix, nativeBlockNum);

	fseek(fp_in, 0L, SEEK_SET);
int fread_size;

//	fread(codeBuf, sizeof(uint8_t), chunkSize, fp_in);
	fread_size = fread(codeBuf, sizeof(uint8_t), chunkSize, fp_in);
	fclose(fp_in);

	fp_in = fopen("native_1", "rb");
//	fseek(fp_in, 0L, SEEK_SET);
	fread_size = fread(codeBuf+chunkSize, sizeof(uint8_t), chunkSize, fp_in);
	fclose(fp_in);

	fp_in = fopen("code_0", "rb");
//	fseek(fp_in, 0L, SEEK_SET);
	fread_size = fread(codeBuf+2*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
	fclose(fp_in);

	fp_in = fopen("code_1", "rb");
//	fseek(fp_in, 0L, SEEK_SET);
	fread_size = fread(codeBuf+3*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
	fclose(fp_in);

	cudaMemcpy(codeBuf_d, codeBuf, codeSize, cudaMemcpyHostToDevice);

/*
int i;
for(i=0; i<4*chunkSize; i++)
{
printf("%c", codeBuf[i]);
}
*/
	int gridDimX = (int)(ceil((float)chunkSize/TILE_WIDTH));
	int gridDimY = (int)(ceil((float)parityBlockNum/TILE_WIDTH));
	dim3 grid(gridDimX, gridDimY);
	dim3 block(TILE_WIDTH, TILE_WIDTH);
	decode_chunk<<<grid, block>>>(dataBuf_d, decodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
	cudaMemcpy(dataBuf, dataBuf_d, dataSize, cudaMemcpyDeviceToHost);

	fp_out = fopen("test", "wb");
	fwrite(dataBuf, sizeof(uint8_t), nativeBlockNum*chunkSize, fp_out);
	fclose(fp_out);

	cudaFree(decodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);

	free(decodingMatrix);
	free(dataBuf);
	free(codeBuf);
#endif
/*
	int matrixSize = nativeBlockNum*nativeBlockNum*sizeof(uint8_t);
	uint8_t *encodingMatrix;	//host
	uint8_t *encodingMatrix_d;	//device
	uint8_t *decodingMatrix;	//host
	uint8_t *decodingMatrix_d;	//device
	encodingMatrix = (uint8_t*) malloc( matrixSize );
	decodingMatrix = (uint8_t*) malloc( matrixSize );
	cudaMalloc( (void **)&encodingMatrix_d, matrixSize );
	cudaMalloc( (void **)&decodingMatrix_d, matrixSize );
	cudaMemcpy(encodingMatrix_d, encodingMatrix, matrixSize, cudaMemcpyHostToDevice);
	GPUGausSeidel (encodingMatrix_d, decodingMatrix_d, nativeBlockNum);
//	GPUGausSeidel (encodingMatrix_d, decodingMatrix_d);
	cudaMemcpy(decodingMatrix, decodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
//	cudaDeviceSynchronize();

	int gridDimX = (int)(ceil((float)chunkSize/TILE_WIDTH));
	int gridDimY = (int)(ceil((float)parityBlockNum/TILE_WIDTH));
	dim3 grid(gridDimX, gridDimY);
	dim3 block(TILE_WIDTH, TILE_WIDTH);
	decode_chunk<<<grid, block>>>(dataBuf_d, decodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
	cudaMemcpy(dataBuf, dataBuf_d, dataSize, cudaMemcpyDeviceToHost);


	char output_file_name[10];
	for(i=0; i<parityBlockNum; i++)
	{
		sprintf(output_file_name, "code_%d", i);
		if( ( fp_out = fopen(output_file_name, "wb") ) == NULL )
		{
			printf("Can not open source file!\n");
			exit(0);
		}
		if( fwrite(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_out ) != sizeof(uint8_t)*chunkSize )
		{
			printf("fwrite error!\n");
			exit(0);
		}
		fclose(fp_out);
	}

	for(i=0; i<nativeBlockNum; i++)
	{
		sprintf(output_file_name, "native_%d", i);
		if( ( fp_out = fopen(output_file_name, "wb") ) == NULL )
		{
			printf("Can not open source file!\n");
			exit(0);
		}
		if( fwrite(dataBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_out ) != sizeof(uint8_t)*chunkSize )
		{
			printf("fwrite error!\n");
			exit(0);
		}
		fclose(fp_out);
	}

//	cudaFree(encodingMatrix_d);
	cudaFree(decodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);

//	free(encodingMatrix);
	free(decodingMatrix);
	free(dataBuf);
	free(codeBuf);
*/

}

