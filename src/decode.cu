/*
 * =====================================================================================
 *
 *       Filename:  decode.cu
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/05/2012 10:50:55 PM
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
#include "matrix.h"

// #define DEBUG

void show_squre_matrix(uint8_t *matrix, int size)
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
}

void copy_matrix(uint8_t *src, uint8_t *des, int srcRowIndex, int desRowIndex, int rowSize)
{
	int i;
	for(i=0; i<rowSize; i++)
	{
		des[desRowIndex*rowSize+i] = src[srcRowIndex*rowSize+i];
	}
}

extern "C"
void decode(uint8_t *dataBuf, uint8_t *codeBuf, uint8_t *encodingMatrix, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	int codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	uint8_t *dataBuf_d;		//device
	uint8_t *codeBuf_d;		//device

	float totalComputationTime = 0;
	float totalCommunicationTime = 0;
	// compute total execution time
	float totalTime;
	cudaEvent_t totalStart, totalStop;
	// create event
	cudaEventCreate(&totalStart);
	cudaEventCreate(&totalStop);
	cudaEventRecord(totalStart);

	cudaMalloc( (void **)&dataBuf_d, dataSize );
//	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc( (void **)&codeBuf_d, codeSize );
//	cudaMemset(codeBuf_d, 0, codeSize);

	// compute step execution time
	float stepTime;
	cudaEvent_t stepStart, stepStop;
	// create event
	cudaEventCreate(&stepStart);
	cudaEventCreate(&stepStop);

	// record event
	cudaEventRecord(stepStart);
	cudaMemcpy(codeBuf_d, codeBuf, codeSize, cudaMemcpyHostToDevice);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("Copy code from CPU to GPU: %fms\n", stepTime);
	totalCommunicationTime += stepTime;

	int matrixSize = nativeBlockNum * nativeBlockNum;
	uint8_t *encodingMatrix_d;	//device
	uint8_t *decodingMatrix_d;	//device
	cudaMalloc( (void **)&encodingMatrix_d, matrixSize );
	cudaMalloc( (void **)&decodingMatrix_d, matrixSize );

	// record event
	cudaEventRecord(stepStart);
	cudaMemcpy(encodingMatrix_d, encodingMatrix, matrixSize, cudaMemcpyHostToDevice);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("Copy encoding matrix from CPU to GPU: %fms\n", stepTime);
	totalCommunicationTime += stepTime;

	// record event
	cudaEventRecord(stepStart);
	invert_matrix(encodingMatrix_d, decodingMatrix_d, nativeBlockNum);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("Generating decoding matrix completed: %fms\n", stepTime);
	totalComputationTime += stepTime;

#ifdef DEBUG
	uint8_t *decodingMatrix;	//host
	decodingMatrix = (uint8_t*) malloc( matrixSize );
	cudaMemcpy(decodingMatrix, decodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
	show_squre_matrix(decodingMatrix, nativeBlockNum);
	free(decodingMatrix);
#endif

//	int gridDimX = (int)(ceil((float)chunkSize/TILE_WIDTH));
//	int gridDimY = (int)(ceil((float)parityBlockNum/TILE_WIDTH));
//	dim3 grid(gridDimX, gridDimY);
//	dim3 block(TILE_WIDTH, TILE_WIDTH);

//	int gridDimX = (int)( ceil((float)chunkSize / TILE_WIDTH_COL) );
	int gridDimX = min( (int)( ceil((float)chunkSize / TILE_WIDTH_COL) ), SINGLE_GRID_SIZE );
	int gridDimY = (int)( ceil((float)nativeBlockNum / TILE_WIDTH_ROW) );
	dim3 grid(gridDimX, gridDimY);
//	dim3 block(TILE_WIDTH_ROW, TILE_WIDTH_COL);
	dim3 block(TILE_WIDTH_COL, TILE_WIDTH_ROW);
	// record event
	cudaEventRecord(stepStart);
	decode_chunk<<<grid, block>>>(dataBuf_d, decodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("Decoding file completed: %fms\n", stepTime);
	totalComputationTime += stepTime;

	// record event
	cudaEventRecord(stepStart);
	cudaMemcpy(dataBuf, dataBuf_d, dataSize, cudaMemcpyDeviceToHost);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("copy data from GPU to CPU: %fms\n", stepTime);
	totalCommunicationTime += stepTime;

	cudaFree(decodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);

	// record event and synchronize
	cudaEventRecord(totalStop);
	cudaEventSynchronize(totalStop);
	// get event elapsed time
	cudaEventElapsedTime(&totalTime, totalStart, totalStop);
	printf("Total computation time: %fms\n", totalComputationTime);
	printf("Total communication time: %fms\n", totalCommunicationTime);
	printf("Total GPU decoding time: %fms\n", totalTime);
}

extern "C"
void decode_file(char *inFile, char *confFile, char *outFile)
{
	int chunkSize = 1;
	int totalSize;
	int parityBlockNum;
	int nativeBlockNum;

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host

	int dataSize;
	int codeSize;

	FILE *fp_meta;
	FILE *fp_in;
	FILE *fp_out;

	int totalMatrixSize;
	int matrixSize;
	uint8_t *totalEncodingMatrix;	//host
	uint8_t *encodingMatrix;	//host
	char metadata_file_name[strlen(inFile) + 15];
	sprintf(metadata_file_name, "%s.METADATA", inFile);
	if((fp_meta = fopen(metadata_file_name, "rb")) == NULL)
	{
		printf("Can not open metadata file!\n");
		exit(0);
	}
	fscanf(fp_meta, "%d", &totalSize);
	fscanf(fp_meta, "%d %d", &parityBlockNum, &nativeBlockNum);

	chunkSize = (totalSize / nativeBlockNum) + ( totalSize%nativeBlockNum != 0 ); 
//	chunkSize = (int) (ceil( (double)totalSize / nativeBlockNum )); 
#ifdef DEBUG
printf("chunk size: %d\n", chunkSize);
#endif

	totalMatrixSize = nativeBlockNum * ( nativeBlockNum + parityBlockNum );
	totalEncodingMatrix = (uint8_t*) malloc( totalMatrixSize );
	matrixSize = nativeBlockNum * nativeBlockNum;
	encodingMatrix = (uint8_t*) malloc( matrixSize );
	for (int i = 0; i < totalMatrixSize; ++i)
	{
		fscanf(fp_meta, "%d", &totalEncodingMatrix[i]);
	}
//	fclose(fp_meta);

	dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( dataSize );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( codeSize);
	memset(codeBuf, 0, codeSize);

	FILE *fp_conf;
	char input_file_name[strlen(inFile) + 20];
	int index;
	fp_conf = fopen(confFile, "r");

	for(int i = 0; i < nativeBlockNum; i++)
	{
		fscanf(fp_conf, "%s", input_file_name);
		index = atoi(input_file_name+1);

		copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

		fp_in = fopen(input_file_name, "rb");
		fseek(fp_in, 0L, SEEK_SET);
		// this part can be process in parallel with computing inversed matrix
		fread(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
		fclose(fp_in);
	}
	fclose(fp_conf);
	
	decode(dataBuf, codeBuf, encodingMatrix, nativeBlockNum, parityBlockNum, chunkSize);

	if(outFile == NULL)
	{
		fp_out = fopen(inFile, "wb");
	}
	else
	{
		fp_out = fopen(outFile, "wb");
	}
	fwrite(dataBuf, sizeof(uint8_t), totalSize, fp_out);
	fclose(fp_out);

	free(dataBuf);
	free(codeBuf);

}
