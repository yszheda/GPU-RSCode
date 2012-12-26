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
#include "matrix.h"

//#define W 8
//#define NW (1 << W) /* In other words, NW equals 2 to the w-th power */

//#define TILE_WIDTH 2
//#define TILE_DEPTH 2

//#define BUFFER_SIZE 256

//#define DEBUG 

//__global__ void encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
//{
//	matrix_mul(parityCoeff, dataChunk, codeChunk, parityBlockNum, nativeBlockNum, chunkSize);
//}
//
//__global__ void gen_encoding_matrix(uint8_t *encodingMatrix, int row, int col)
//{
//	int i = threadIdx.x;
//	int j = threadIdx.y;
//	setup_tables(8);
//	__syncthreads();
//	encodingMatrix[i*col + j] = gf_pow(j+1, i);
//}


void write_metadata(int totalSize, int parityBlockNum, int nativeBlockNum, uint8_t *encodingMatrix)
{
	FILE *fp;
	if( ( fp = fopen(".METADATA", "wb") ) == NULL )
	{
		printf("Can not open META file!\n");
		exit(0);
	}
	fprintf(fp, "%d\n", totalSize);
	fprintf(fp, "%d %d\n", parityBlockNum, nativeBlockNum);
	int i;
	int j;
	for(i=0; i<nativeBlockNum; i++)
	{
		for(j=0; j<nativeBlockNum; j++)
		{
			if(i == j)
			{
				fprintf(fp, "1 ");
			}
			else
			{
				fprintf(fp, "0 ");
			}
		}
		fprintf(fp, "\n");
	}
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

extern "C"
void encode(uint8_t *dataBuf, uint8_t *codeBuf, int nativeBlockNum, int parityBlockNum, int chunkSize, int totalSize)
{
	uint8_t *dataBuf_d;		//device
	uint8_t *codeBuf_d;		//device
	int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	int codeSize = parityBlockNum*chunkSize*sizeof(uint8_t);
	cudaMalloc( (void **)&dataBuf_d, nativeBlockNum*chunkSize*sizeof(uint8_t) );
//	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc( (void **)&codeBuf_d, parityBlockNum*chunkSize*sizeof(uint8_t) );
//	cudaMemset(codeBuf_d, 0, codeSize);

	cudaMemcpy(dataBuf_d, dataBuf, nativeBlockNum*chunkSize*sizeof(uint8_t), cudaMemcpyHostToDevice);

	uint8_t *encodingMatrix;	//host
	uint8_t *encodingMatrix_d;	//device
	encodingMatrix = (uint8_t*) malloc( parityBlockNum*nativeBlockNum*sizeof(uint8_t) );
	cudaMalloc( (void **)&encodingMatrix_d, parityBlockNum*nativeBlockNum*sizeof(uint8_t) );

	dim3 blk(parityBlockNum, nativeBlockNum);
	gen_encoding_matrix<<<1, blk>>>(encodingMatrix_d, parityBlockNum, nativeBlockNum);
//	cudaDeviceSynchronize();

	cudaMemcpy(encodingMatrix, encodingMatrix_d, parityBlockNum*nativeBlockNum*sizeof(uint8_t), cudaMemcpyDeviceToHost);
	write_metadata(totalSize, parityBlockNum, nativeBlockNum, encodingMatrix);

//	int gridDimX = (int)(ceil((float)chunkSize/TILE_WIDTH));
//	int gridDimY = (int)(ceil((float)parityBlockNum/TILE_WIDTH));
//	dim3 grid(gridDimX, gridDimY);
//	dim3 block(TILE_WIDTH, TILE_WIDTH);
	int gridDimX = (int)( ceil((float)chunkSize / TILE_WIDTH_COL) );
	int gridDimY = (int)( ceil((float)nativeBlockNum / TILE_WIDTH_ROW) );
	dim3 grid(gridDimX, gridDimY);
	dim3 block(TILE_WIDTH_ROW, TILE_WIDTH_COL);
	encode_chunk<<<grid, block>>>(dataBuf_d, encodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);

	cudaMemcpy(codeBuf, codeBuf_d, parityBlockNum*chunkSize*sizeof(uint8_t), cudaMemcpyDeviceToHost);

	cudaFree(encodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);

	free(encodingMatrix);

}

extern "C"
void encode_file(char *file, int nativeBlockNum, int parityBlockNum)
{
	int chunkSize = 1;
	int totalSize;

	FILE *fp_in;
	FILE *fp_out;
	if( ( fp_in = fopen(file,"rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}

	fseek(fp_in, 0L, SEEK_END);
	//ftell() get the total size of the file
	totalSize = ftell(fp_in);
	chunkSize = (int) (ceil( (float)ftell(fp_in) / nativeBlockNum )); 

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host
	int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	int codeSize = parityBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( nativeBlockNum*chunkSize*sizeof(uint8_t) );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( parityBlockNum*chunkSize*sizeof(uint8_t) );
	memset(codeBuf, 0, codeSize);
	
	int i;
	for(i=0; i<nativeBlockNum; i++)
	{
		if( fseek( fp_in, i*chunkSize, SEEK_SET ) == -1 )
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
	
	encode(dataBuf, codeBuf, nativeBlockNum, parityBlockNum, chunkSize, totalSize);

	char output_file_name[20];
	for(i=0; i<nativeBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_", i);
		strcat(output_file_name, file);
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
	for(i=0; i<parityBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_", i+nativeBlockNum);
		strcat(output_file_name, file);
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

	free(dataBuf);
	free(codeBuf);
}
/*
int main(int argc, char *argv[])
{
	int nativeBlockNum = 4;
	int parityBlockNum = 2;
//	int chunkSize = sizeof(uint8_t);
	int chunkSize = 1;
	int totalSize;

	FILE *fp_in;
	FILE *fp_out;
	if( ( fp_in = fopen(argv[1],"rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}

	fseek(fp_in, 0L, SEEK_END);
	//ftell() get the total size of the file
	totalSize = ftell(fp_in);
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
	dim3 blk(parityBlockNum, nativeBlockNum);
	gen_encoding_matrix<<<1, blk>>>(encodingMatrix_d, parityBlockNum, nativeBlockNum);
//	cudaDeviceSynchronize();

	cudaMemcpy(encodingMatrix, encodingMatrix_d, parityBlockNum*nativeBlockNum*sizeof(uint8_t), cudaMemcpyDeviceToHost);
	write_metadata(totalSize, parityBlockNum, nativeBlockNum, encodingMatrix);

	int gridDimX = (int)(ceil((float)chunkSize/TILE_WIDTH));
	int gridDimY = (int)(ceil((float)parityBlockNum/TILE_WIDTH));
	dim3 grid(gridDimX, gridDimY);
	dim3 block(TILE_WIDTH, TILE_WIDTH);

	encode_chunk<<<grid, block>>>(dataBuf_d, encodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
	cudaMemcpy(codeBuf, codeBuf_d, parityBlockNum*chunkSize*sizeof(uint8_t), cudaMemcpyDeviceToHost);


#ifdef DEBUG
//	show_code_chunk(codeBuf, parityBlockNum, chunkSize);
#endif

	char output_file_name[20];

	for(i=0; i<nativeBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_", i);
		strcat(output_file_name, argv[1]);
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

	for(i=0; i<parityBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_", i+nativeBlockNum);
		strcat(output_file_name, argv[1]);
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

	cudaFree(encodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);

	free(encodingMatrix);
	free(dataBuf);
	free(codeBuf);

}
*/
