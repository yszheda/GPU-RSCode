/*
 * =====================================================================================
 *
 *       Filename:  test-seq.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/05/2012 10:42:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
//#include "galoisfield.h"

//extern uint8_t encodingMatrix[];
//extern unsigned char encodingMatrix[];
unsigned char* encodingMatrix;

const int w = 8;

const int NW = 1 << 8;
//#define NW (1 << w) /* In other words, NW equals 2 to the w-th power */
//#define DEBUG 
#define BUFFER_SIZE 256

//unsigned int prim_poly_4 = 023;
//unsigned int prim_poly_8 = 0435;
//unsigned int prim_poly_16 = 0210013;
unsigned int prim_poly_4 = 023;
unsigned int prim_poly_8 = 0435;
//uint8_t prim_poly_8 = 285;
unsigned int prim_poly_16 = 0210013;
uint8_t *gflog;
uint8_t *gfilog;

int setup_tables(int w)
{
	unsigned int b;
	unsigned int r;
   	unsigned int log;
	unsigned int x_to_w;
	unsigned int prim_poly;
	unsigned int x;
	unsigned int y;
//	uint8_t b, log, x_to_w, prim_poly;
	switch(w) 
	{
		case 4: prim_poly = prim_poly_4; break;
		case 8: prim_poly = prim_poly_8; break;
		case 16: prim_poly = prim_poly_16; break;
		default: return -1;
	}
	x_to_w = 1 << w;
	gflog = (uint8_t *) malloc (sizeof(uint8_t) * x_to_w);
	gfilog = (uint8_t *) malloc (sizeof(uint8_t) * x_to_w);
	b = 1;
	r = 0;
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
		gfilog[log] = (uint8_t) b;
		b = b << 1;
		if (b & x_to_w) 
		{
			b = b ^ prim_poly;
//#ifdef DEBUG
//printf("prim_poly=%d\n", prim_poly);
//printf("test b=%d\n", b);
//#endif
		}
	}
//#ifdef DEBUG
//	printf("b=%d\n",b);
//#endif
	return 0;
}

void showgflog()
{
	int i;
	printf("gflog\n");
	for(i =0; i< NW; i++)
	{
		printf("%d ", gflog[i]);
	}
}

void showgfilog()
{
	int i;
	printf("gfilog\n");
	for(i =0; i< NW; i++)
	{
		printf("%d ", gfilog[i]);
	}
}

uint8_t gf_add(uint8_t a, uint8_t b)
{
	return a^b;
}

uint8_t gf_sub(uint8_t a, uint8_t b)
{
	return gf_add(a, b);
}

uint8_t gf_mul(uint8_t a, uint8_t b)
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
	return gfilog[sum_log];
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

uint8_t gf_div(uint8_t a, uint8_t b)
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
	return gfilog[diff_log];
}

uint8_t gf_pow(uint8_t a, uint8_t power)
{
	uint8_t pow_log = (gflog[a] * power) % (NW-1);
	return gfilog[pow_log];
}

void gen_encoding_matrix(int row, int col)
{
	int i;
	int j;
	encodingMatrix = (unsigned char*) malloc( row*col*sizeof(unsigned char) );
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			encodingMatrix[i*col + j] = gf_pow(j+1, i);
		}
	}
}

// C=AB
// A: nxp
// B: pxm
// C: nxm
void matrix_mul(unsigned char *A, unsigned char *B, unsigned char *C, int n, int p, int m)
{
	int i;
	int j;
	int k;
	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
		{
			for(k=0; k<p; k++)
			{
				C[i*m+j] = gf_add(C[i*m+j], gf_mul(A[i*p+k],B[k*m+j]));
			}
		}
	}
}

void encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
//	codeChunk = (unsigned char*)malloc(parityBlockNum*chunkSize);
	matrix_mul(parityCoeff, dataChunk, codeChunk, parityBlockNum, nativeBlockNum, chunkSize);
}

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

int main(int argc, char *argv[])
{
	int nativeBlockNum = 4;
	int parityBlockNum = 2;
//	int chunkSize = sizeof(uint8_t);
	int chunkSize = 1;

//	char[] file_name_src = argv[1];
	FILE *fp_in;
	FILE *fp_out;
	int i;
	if( ( fp_in = fopen(argv[1],"rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}

	fseek(fp_in, 0L, SEEK_END);
	chunkSize = ftell(fp_in)/nativeBlockNum; //ftell() get the total size of the file

	uint8_t *dataBuf;
	uint8_t *codeBuf;
	dataBuf = (unsigned char*) malloc( nativeBlockNum*chunkSize*sizeof(uint8_t) );
	codeBuf = (unsigned char*) malloc( parityBlockNum*chunkSize*sizeof(uint8_t) );

	for(i=0; i<nativeBlockNum; i++)
	{
		if( fseek(fp_in, i*chunkSize, SEEK_SET) == -1 )
		{
			printf("fseek error!\n");
			exit(0);
		}

		if( fread(dataBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in ) != sizeof(uint8_t)*chunkSize )
		{
			printf("fread error!\n");
			exit(0);
		}
	}
	fclose(fp_in);

	// setup table for GF(2^8)
	setup_tables(8);
	gen_encoding_matrix(parityBlockNum, nativeBlockNum);
	encode_chunk(dataBuf, encodingMatrix, codeBuf, nativeBlockNum, parityBlockNum, chunkSize);

#ifdef DEBUG
	show_code_chunk(codeBuf, parityBlockNum, chunkSize);
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



//	uint8_t testData[4] = {9,1,2,0};
//	uint8_t *testCode;
//
//	// setup table for GF(2^8)
//	setup_tables(8);
//	gen_encoding_matrix(parityBlockNum, nativeBlockNum);
//	testCode = (unsigned char*)malloc(parityBlockNum*chunkSize);
//	encode_chunk(testData, encodingMatrix, testCode, nativeBlockNum, parityBlockNum, chunkSize);
//
//#ifdef DEBUG
//showgflog();
//printf("\n");
//showgfilog();
//printf("\n");
//#endif
//	show_code_chunk(testCode, parityBlockNum, chunkSize);
//#ifdef DEBUG
//printf("\n");
//show_encoding_matrix(parityBlockNum, nativeBlockNum);
//#endif

free(encodingMatrix);

}
