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


	float time;
	// compute the execution time
	cudaEvent_t start, stop;
	// create event
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	// record event
	cudaEventRecord(start);
	dim3 blk(parityBlockNum, nativeBlockNum);
	gen_encoding_matrix<<<1, blk>>>(encodingMatrix_d, parityBlockNum, nativeBlockNum);
//	cudaDeviceSynchronize();
	// record event and synchronize
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	// get event elapsed time
	cudaEventElapsedTime(&time, start, stop);
	printf("Generating encoding matrix completed: %fms\n", time);

	cudaMemcpy(encodingMatrix, encodingMatrix_d, parityBlockNum*nativeBlockNum*sizeof(uint8_t), cudaMemcpyDeviceToHost);
	write_metadata(totalSize, parityBlockNum, nativeBlockNum, encodingMatrix);

	// record event
	cudaEventRecord(start);
//	int gridDimX = (int)(ceil((float)chunkSize/TILE_WIDTH));
//	int gridDimY = (int)(ceil((float)parityBlockNum/TILE_WIDTH));
//	dim3 grid(gridDimX, gridDimY);
//	dim3 block(TILE_WIDTH, TILE_WIDTH);
	int gridDimX = min( (int)( ceil((float)chunkSize / TILE_WIDTH_COL) ), SINGLE_GRID_SIZE );
	int gridDimY = (int)( ceil((float)nativeBlockNum / TILE_WIDTH_ROW) );
	dim3 grid(gridDimX, gridDimY);
//	dim3 block(TILE_WIDTH_ROW, TILE_WIDTH_COL);
	dim3 block(TILE_WIDTH_COL, TILE_WIDTH_ROW);
	encode_chunk<<<grid, block>>>(dataBuf_d, encodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
	// record event and synchronize
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	// get event elapsed time
	cudaEventElapsedTime(&time, start, stop);
	printf("Encoding file completed: %fms\n", time);

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
	chunkSize = (ftell(fp_in) / nativeBlockNum) + ( ftell(fp_in)%nativeBlockNum != 0 ); 
//	chunkSize = (int) (ceil( (long double) (ftell(fp_in) / nativeBlockNum)) ); 

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

	char output_file_name[100];
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
