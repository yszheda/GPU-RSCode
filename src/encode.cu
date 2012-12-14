/*
 * =====================================================================================
 *
 *       Filename:  encode.cu
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/05/2012 10:42:32 PM
 *       Revision:  none
 *       Compiler:  nvcc/ggc
 *
 *         Author:  Shuai YUAN (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <stdint.h>
//#include "galoisfield.h"

//const int w = 8;
//const int NW = 1 << 8;

#define W 8
#define NW (1 << W) /* In other words, NW equals 2 to the w-th power */

//#define DEBUG 
#define BUFFER_SIZE 256

__shared__ uint8_t gflog[512];
__shared__ uint8_t gfexp[512];

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

__global__ void gen_encoding_matrix(uint8_t *encodingMatrix, int row, int col)
{
	int i = threadIdx.x;
	int j = threadIdx.y;
	setup_tables(8);
	encodingMatrix[i*col + j] = gf_pow(j+1, i);
}

// C=AB
// A: nxp
// B: pxm
// C: nxm
__device__ void matrix_mul(unsigned char *A, unsigned char *B, unsigned char *C, int n, int p, int m)
{
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
}

__global__ void encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	matrix_mul(parityCoeff, dataChunk, codeChunk, parityBlockNum, nativeBlockNum, chunkSize);
}
/*
void show_code_chunk(unsigned char *codeChunk, int parityBlockNum, int chunkSize)
{
	int i;
	int j;
	for(i=0; i<parityBlockNum; i++)
	{
		for(j=0; j<chunkSize; j++)
		{
			printf("%d ", codeChunk[i*chunkSize+j]);
		}
		printf("\n");
	}
}

void show_encoding_matrix(int parityBlockNum, int nativeBlockNum)
{
	int i;
	int j;
	for(i=0; i<parityBlockNum; i++)
	{
		for(j=0; j<nativeBlockNum; j++)
		{
			printf("%d ", encodingMatrix[i*nativeBlockNum+j]);
		}
		printf("\n");
	}
}
*/

void write_metadata(int parityBlockNum, int nativeBlockNum, uint8_t *encodingMatrix)
{
	FILE *fp;
	if( ( fp = fopen(".METADATA", "wb") ) == NULL )
	{
		printf("Can not open META file!\n");
		exit(0);
	}
	fprintf(fp, "%d %d\n", parityBlockNum, nativeBlockNum);
	int i;
	int j;
	for(i=0; i<parityBlockNum; i++)
	{
		for(j=0; j<nativeBlockNum; j++)
		{
			fprintf(fp, "%d ", encodingMatrix[i*nativeBlockNum+j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

int main(int argc, char *argv[])
{
	int nativeBlockNum = 4;
	int parityBlockNum = 2;
//	int chunkSize = sizeof(uint8_t);
	int chunkSize = 1;

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

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host
	uint8_t *dataBuf_d;		//device
	uint8_t *codeBuf_d;		//device
	int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	int codeSize = parityBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( nativeBlockNum*chunkSize*sizeof(uint8_t) );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( parityBlockNum*chunkSize*sizeof(uint8_t) );
	memset(codeBuf, 0, codeSize);
	cudaMalloc( (void **)&dataBuf_d, nativeBlockNum*chunkSize*sizeof(uint8_t) );
	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc( (void **)&codeBuf_d, parityBlockNum*chunkSize*sizeof(uint8_t) );
	cudaMemset(codeBuf_d, 0, codeSize);

	int i;
	for(i=0; i<nativeBlockNum; i++)
	{
		if( fseek( fp_in, i*chunkSize, SEEK_SET) == -1 )
		{
			printf("fseek error!\n");
			exit(0);
		}

		if( fread( dataBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in ) == EOF )
		{
			printf("fread error!\n");
			exit(0);
		}
	}
	fclose(fp_in);

	cudaMemcpy(dataBuf_d, dataBuf, nativeBlockNum*chunkSize*sizeof(uint8_t), cudaMemcpyHostToDevice);

	// setup table for GF(2^8)
//	uint8_t *gflog;	//host
//	uint8_t *gfexp;	//host
//	uint8_t *gflog_d;	//device
//	uint8_t *gfexp_d;	//device
//	gflog = (uint8_t*) malloc( 512*sizeof(uint8_t) );
//	gfexp = (uint8_t*) malloc( 512*sizeof(uint8_t) );
//	cudaMalloc( (void **)&gflog_d, 512*sizeof(uint8_t) );
//	cudaMalloc( (void **)&gfexp_d, 512*sizeof(uint8_t) );
//	setup_tables(8);
//	cudaDeviceSynchronize();
//	cudaMemcpy(gflog_d, gflog, 512*sizeof(uint8_t), cudaMemcpyDeviceToHost);

	uint8_t *encodingMatrix;	//host
	uint8_t *encodingMatrix_d;	//device
	encodingMatrix = (uint8_t*) malloc( parityBlockNum*nativeBlockNum*sizeof(uint8_t) );
	cudaMalloc( (void **)&encodingMatrix_d, parityBlockNum*nativeBlockNum*sizeof(uint8_t) );
	dim3 block(parityBlockNum, nativeBlockNum);
	gen_encoding_matrix<<<1, block>>>(encodingMatrix_d, parityBlockNum, nativeBlockNum);
//	cudaDeviceSynchronize();

	cudaMemcpy(encodingMatrix, encodingMatrix_d, parityBlockNum*nativeBlockNum*sizeof(uint8_t), cudaMemcpyDeviceToHost);
	write_metadata(parityBlockNum, nativeBlockNum, encodingMatrix);

	encode_chunk<<<1, 1>>>(dataBuf_d, encodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
	cudaMemcpy(codeBuf, codeBuf_d, parityBlockNum*chunkSize*sizeof(uint8_t), cudaMemcpyDeviceToHost);


#ifdef DEBUG
//	show_code_chunk(codeBuf, parityBlockNum, chunkSize);
#endif

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

	cudaFree(encodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);

	free(encodingMatrix);
	free(dataBuf);
	free(codeBuf);

}
