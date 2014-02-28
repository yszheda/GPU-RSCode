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
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "matrix.h"

void write_metadata(char *fileName, int totalSize, int parityBlockNum, int nativeBlockNum, uint8_t* encodingMatrix)
{
	FILE *fp;
	if( ( fp = fopen(fileName, "wb") ) == NULL )
	{
		printf("Can not open META file!\n");
		exit(0);
	}
	fprintf(fp, "%d\n", totalSize);
	fprintf(fp, "%d %d\n", parityBlockNum, nativeBlockNum);
	for (int i = 0; i < nativeBlockNum; ++i)
	{
		for (int j = 0; j < nativeBlockNum; ++j)
		{
			if (i == j)
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
	for (int i = 0; i < parityBlockNum; ++i)
	{
		for (int j = 0; j < nativeBlockNum; ++j)
		{
			fprintf(fp, "%d ", encodingMatrix[i*nativeBlockNum + j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

extern "C"
void encode(char *fileName, uint8_t *dataBuf, uint8_t *codeBuf, int nativeBlockNum, int parityBlockNum, int chunkSize, int totalSize, int gridDimXSize, int streamNum)
{
	int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	int codeSize = parityBlockNum*chunkSize*sizeof(uint8_t);

	float totalComputationTime = 0;
	float totalCommunicationTime = 0;
	// compute total execution time
	float totalTime;
	cudaEvent_t totalStart, totalStop;
	// create event
	cudaEventCreate(&totalStart);
	cudaEventCreate(&totalStop);
	cudaEventRecord(totalStart);

	// compute step execution time
	float stepTime;
	cudaEvent_t stepStart, stepStop;
	// create event
	cudaEventCreate(&stepStart);
	cudaEventCreate(&stepStop);

	// record event
//	cudaEventRecord(stepStart);
//	cudaMemcpy(dataBuf_d, dataBuf, nativeBlockNum*chunkSize*sizeof(uint8_t), cudaMemcpyHostToDevice);
//	// record event and synchronize
//	cudaEventRecord(stepStop);
//	cudaEventSynchronize(stepStop);
//	// get event elapsed time
//	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
//	printf("Copy data from CPU to GPU: %fms\n", stepTime);
//	totalCommunicationTime += stepTime;

	uint8_t *encodingMatrix;	//host
	uint8_t *encodingMatrix_d;	//device
//	encodingMatrix = (uint8_t*) malloc( parityBlockNum*nativeBlockNum*sizeof(uint8_t) );
	cudaMallocHost( (void **)&encodingMatrix, parityBlockNum*nativeBlockNum*sizeof(uint8_t) );
	cudaMalloc( (void **)&encodingMatrix_d, parityBlockNum*nativeBlockNum*sizeof(uint8_t) );

	// record event
	cudaEventRecord(stepStart);
	dim3 blk(parityBlockNum, nativeBlockNum);
	gen_encoding_matrix<<<1, blk>>>(encodingMatrix_d, parityBlockNum, nativeBlockNum);
//	cudaDeviceSynchronize();
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("Generating encoding matrix completed: %fms\n", stepTime);
	totalComputationTime += stepTime;

	// record event
	cudaEventRecord(stepStart);
	cudaMemcpy(encodingMatrix, encodingMatrix_d, parityBlockNum*nativeBlockNum*sizeof(uint8_t), cudaMemcpyDeviceToHost);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("Copy encoding matrix from GPU to CPU: %fms\n", stepTime);
	totalCommunicationTime += stepTime;

    // use cuda stream to encode the file
    // to obtain computation and comunication overhead
//    int streamNum = ( chunkSize / STREAM_MAX_CHUNK_SIZE ) + ( chunkSize%STREAM_MAX_CHUNK_SIZE != 0 ); 
//    int streamChunkSize = min(chunkSize, STREAM_MAX_CHUNK_SIZE);

	int streamMaxChunkSize = (chunkSize / streamNum) + (chunkSize % streamNum != 0);
    cudaStream_t stream[streamNum];
    //NOTE: need cudaMallocHost
    for(int i = 0; i < streamNum; i++)
    {
		cudaStreamCreate(&stream[i]);
    }

	uint8_t *dataBuf_d[streamNum];		//device
	uint8_t *codeBuf_d[streamNum];		//device
    for(int i = 0; i < streamNum; i++)
    {
        int streamChunkSize = min(chunkSize-i*streamMaxChunkSize, streamMaxChunkSize);

//        uint8_t *dataBuf_d;                //device
//        uint8_t *codeBuf_d;                //device
        int dataSize = nativeBlockNum*streamChunkSize*sizeof(uint8_t);
        int codeSize = parityBlockNum*streamChunkSize*sizeof(uint8_t);

        cudaMalloc( (void **)&dataBuf_d[i], dataSize );
        cudaMalloc( (void **)&codeBuf_d[i], codeSize );
	}

    for(int i = 0; i < streamNum; i++)
    {
        int streamChunkSize = min(chunkSize-i*streamMaxChunkSize, streamMaxChunkSize);
        int dataSize = nativeBlockNum*streamChunkSize*sizeof(uint8_t);
        int codeSize = parityBlockNum*streamChunkSize*sizeof(uint8_t);
//        // record event
//        cudaEventRecord(stepStart);
        for(int j = 0; j < nativeBlockNum; j++)
        {
			cudaMemcpyAsync(dataBuf_d[i]+j*streamChunkSize, dataBuf+j*chunkSize+i*streamChunkSize, 
                                                streamChunkSize*sizeof(uint8_t), 
                                                cudaMemcpyHostToDevice, stream[i]);
        }
//        // record event and synchronize
//        cudaEventRecord(stepStop);
//        cudaEventSynchronize(stepStop);
//        // get event elapsed time
//        cudaEventElapsedTime(&stepTime, stepStart, stepStop);
//        printf("Copy data from CPU to GPU in stream: %fms\n", stepTime);
//        totalCommunicationTime += stepTime;

		stepTime = encode_chunk(dataBuf_d[i], encodingMatrix_d, codeBuf_d[i], nativeBlockNum, parityBlockNum, streamChunkSize, stream[i]);
//        printf("Encoding file in stream completed: %fms\n", stepTime);
//        totalComputationTime += stepTime;
		
//        // record event
//        cudaEventRecord(stepStart);
        for(int j = 0; j < parityBlockNum; j++)
        {
                cudaMemcpyAsync(codeBuf+j*chunkSize+i*streamChunkSize, codeBuf_d[i]+j*streamChunkSize, 
                                                streamChunkSize*sizeof(uint8_t),
                                                cudaMemcpyDeviceToHost, stream[i]);
        }
//        // record event and synchronize
//        cudaEventRecord(stepStop);
//        cudaEventSynchronize(stepStop);
//        // get event elapsed time
//        cudaEventElapsedTime(&stepTime, stepStart, stepStop);
//        printf("copy code from GPU to CPU in stream: %fms\n", stepTime);
//        totalCommunicationTime += stepTime;

	}

    for(int i = 0; i < streamNum; i++)
	{
        cudaFree(dataBuf_d[i]);
        cudaFree(codeBuf_d[i]);
    }
/*
    for(int i = 0; i < streamNum; i++)
    {
		cudaStreamDestroy(stream[i]);
    }
*/
	cudaFree(encodingMatrix_d);

	// record event and synchronize
	cudaEventRecord(totalStop);
	cudaEventSynchronize(totalStop);
	// get event elapsed time
	cudaEventElapsedTime(&totalTime, totalStart, totalStop);
//	printf("Total computation time: %fms\n", totalComputationTime);
//	printf("Total communication time: %fms\n", totalCommunicationTime);
	printf("Total GPU encoding time: %fms\n", totalTime);

    for(int i = 0; i < streamNum; i++)
    {
		cudaStreamDestroy(stream[i]);
    }

	char metadata_file_name[strlen(fileName) + 15];
	sprintf(metadata_file_name, "%s.METADATA", fileName);
	write_metadata(metadata_file_name, totalSize, parityBlockNum, nativeBlockNum, encodingMatrix);
//	free(encodingMatrix);
	cudaFreeHost(encodingMatrix);
}

extern "C"
void encode_file(char *fileName, int nativeBlockNum, int parityBlockNum, int gridDimXSize, int streamNum)
{
	int chunkSize = 1;
	int totalSize;

	FILE *fp_in;
	FILE *fp_out;
	if( ( fp_in = fopen(fileName,"rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}

	fseek(fp_in, 0L, SEEK_END);
	// ftell() get the total size of the file
	totalSize = ftell(fp_in);
	chunkSize = (totalSize / nativeBlockNum) + ( totalSize%nativeBlockNum != 0 ); 
//	chunkSize = (ftell(fp_in) / nativeBlockNum) + ( ftell(fp_in)%nativeBlockNum != 0 ); 
//	chunkSize = (int) (ceil( (long double) (ftell(fp_in) / nativeBlockNum)) ); 

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host
	int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	int codeSize = parityBlockNum*chunkSize*sizeof(uint8_t);

//	dataBuf = (uint8_t*) malloc( nativeBlockNum*chunkSize*sizeof(uint8_t) );
//	memset(dataBuf, 0, dataSize);
//	codeBuf = (uint8_t*) malloc( parityBlockNum*chunkSize*sizeof(uint8_t) );
//	memset(codeBuf, 0, codeSize);

	// cuda stream may require cudaMallocHost rather than malloc
	cudaMallocHost( (void **)&dataBuf, dataSize);
	cudaMallocHost( (void **)&codeBuf, codeSize);
	
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
	
	cudaSetDevice(0);
	cudaDeviceProp deviceProperties;
	cudaGetDeviceProperties(&deviceProperties, 0);
	int maxGridDimXSize = deviceProperties.maxGridSize[0];
	if (gridDimXSize > maxGridDimXSize)
	{
		printf("max X dimension grid size is only %d!\n", maxGridDimXSize);
		gridDimXSize = maxGridDimXSize;
	}
	if (gridDimXSize <= 0)
	{
		gridDimXSize = maxGridDimXSize;
	}
	
//	cudaSetDevice(1);

	encode(fileName, dataBuf, codeBuf, nativeBlockNum, parityBlockNum, chunkSize, totalSize, gridDimXSize, streamNum);

	char output_file_name[strlen(fileName) + 5];
	for(i=0; i<nativeBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_%s", i, fileName);
		if( ( fp_out = fopen(output_file_name, "wb") ) == NULL )
		{
			printf("Can not open output file!\n");
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
		sprintf(output_file_name, "_%d_%s", i + nativeBlockNum, fileName);
		if( ( fp_out = fopen(output_file_name, "wb") ) == NULL )
		{
			printf("Can not open output file!\n");
			exit(0);
		}
		if( fwrite(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_out ) != sizeof(uint8_t)*chunkSize )
		{
			printf("fwrite error!\n");
			exit(0);
		}
		fclose(fp_out);
	}

	cudaDeviceReset();

	cudaFreeHost(dataBuf);
	cudaFreeHost(codeBuf);
}
