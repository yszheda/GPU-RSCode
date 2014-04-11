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
#include <pthread.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "matrix.h"

struct ThreadDataType {
	int id;
	int nativeBlockNum;
	int parityBlockNum;
	int chunkSize;
	int totalSize;
	int gridDimXSize;
	int streamNum;
	char* fileName;
	uint8_t* dataBuf;
	uint8_t* codeBuf;
};	/* ----------  end of struct ThreadDataType  ---------- */

typedef struct ThreadDataType ThreadDataType;

static pthread_barrier_t barrier;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  write_metadata
 *  Description:  write metadata into the file with <filename>
 *  The metadata contains:
 *  totalSize: total file size
 *  parityBlockNum: number of parity code chunks
 *  nativeBlockNum: number of native data chunks
 *  encodingMatrix: encoding matrix
 *  The identity matrix will be appended to the encoding matrix as well.
 * =====================================================================================
 */
void write_metadata(char *fileName, int totalSize, int parityBlockNum, int nativeBlockNum, uint8_t* encodingMatrix)
{
	FILE *fp;
	if ((fp = fopen(fileName, "wb")) == NULL)
	{
		printf("Cannot open META file!\n");
		exit(0);
	}
	fprintf(fp, "%d\n", totalSize);
	fprintf(fp, "%d %d\n", parityBlockNum, nativeBlockNum);
	// write the identity matrix into metadata file
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
	// write the encoding matrix into metadata file
	for (int i = 0; i < parityBlockNum; ++i)
	{
		for (int j = 0; j < nativeBlockNum; ++j)
		{
			fprintf(fp, "%d ", encodingMatrix[i * nativeBlockNum + j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  encode
 *  Description:  encode the given buffer of data chunks in the GPU with <id>
 * =====================================================================================
 */
extern "C"
void encode(char *fileName, uint8_t *dataBuf, uint8_t *codeBuf, int id, int nativeBlockNum, int parityBlockNum, int chunkSize, int totalSize, int gridDimXSize, int streamNum)
{
//	cudaSetDevice(id);

	int dataSize = nativeBlockNum * chunkSize * sizeof(uint8_t);
	int codeSize = parityBlockNum * chunkSize * sizeof(uint8_t);

	float totalComputationTime = 0;
	float totalCommunicationTime = 0;
	// compute total execution time
	float totalTime;
	cudaEvent_t totalStart, totalStop;
	// create CUDA event
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
	int matrixSize = parityBlockNum * nativeBlockNum * sizeof(uint8_t);
	// Pageable Host Memory is preferred here since the encodingMatrix is small
	encodingMatrix = (uint8_t*) malloc(matrixSize);
//	cudaMallocHost((void **)&encodingMatrix, matrixSize);
	cudaMalloc((void **)&encodingMatrix_d, matrixSize);

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
	cudaMemcpy(encodingMatrix, encodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("Copy encoding matrix from GPU to CPU: %fms\n", stepTime);
	totalCommunicationTime += stepTime;

	// Use cuda stream to encode the file
	// to achieve computation and comunication overlapping
	// Use DFS way
//	int streamMaxChunkSize = (chunkSize / streamNum) + (chunkSize % streamNum != 0);
	int streamMinChunkSize = chunkSize / streamNum;
	cudaStream_t stream[streamNum];
	for (int i = 0; i < streamNum; i++)
	{
		cudaStreamCreate(&stream[i]);
	}

	uint8_t *dataBuf_d[streamNum];		//device
	uint8_t *codeBuf_d[streamNum];		//device
	for (int i = 0; i < streamNum; i++)
	{
		int streamChunkSize = streamMinChunkSize;
		if (i == streamNum - 1) 
		{
			streamChunkSize = chunkSize - i * streamMinChunkSize;
		} 

		int dataSize = nativeBlockNum * streamChunkSize * sizeof(uint8_t);
		int codeSize = parityBlockNum * streamChunkSize * sizeof(uint8_t);
		
		cudaMalloc((void **)&dataBuf_d[i], dataSize);
		cudaMalloc((void **)&codeBuf_d[i], codeSize);
	}

	for (int i = 0; i < streamNum; i++)
	{
		int streamChunkSize = streamMinChunkSize;
		if (i == streamNum - 1) 
		{
			streamChunkSize = chunkSize - i * streamMinChunkSize;
		} 
		int dataSize = nativeBlockNum * streamChunkSize * sizeof(uint8_t);
		int codeSize = parityBlockNum * streamChunkSize * sizeof(uint8_t);
//		// record event
//		cudaEventRecord(stepStart);
		for (int j = 0; j < nativeBlockNum; j++)
		{
			cudaMemcpyAsync(dataBuf_d[i] + j * streamChunkSize, 
							dataBuf + j * chunkSize + i * streamMinChunkSize, 
							streamChunkSize * sizeof(uint8_t), 
							cudaMemcpyHostToDevice, 
							stream[i]);
		}
//		// record event and synchronize
//		cudaEventRecord(stepStop);
//		cudaEventSynchronize(stepStop);
//		// get event elapsed time
//		cudaEventElapsedTime(&stepTime, stepStart, stepStop);
//		printf("Copy data from CPU to GPU in stream: %fms\n", stepTime);
//		totalCommunicationTime += stepTime;

		stepTime = encode_chunk(dataBuf_d[i], encodingMatrix_d, codeBuf_d[i], nativeBlockNum, parityBlockNum, streamChunkSize, gridDimXSize, stream[i]);
//		printf("Encoding file in stream completed: %fms\n", stepTime);
//		totalComputationTime += stepTime;
		
//		// record event
//		cudaEventRecord(stepStart);
		for (int j = 0; j < parityBlockNum; j++)
		{
			cudaMemcpyAsync(codeBuf + j * chunkSize + i * streamMinChunkSize, 
							codeBuf_d[i] + j * streamChunkSize, 
							streamChunkSize * sizeof(uint8_t),
							cudaMemcpyDeviceToHost, 
							stream[i]);
		}
//		// record event and synchronize
//		cudaEventRecord(stepStop);
//		cudaEventSynchronize(stepStop);
//		// get event elapsed time
//		cudaEventElapsedTime(&stepTime, stepStart, stepStop);
//		printf("copy code from GPU to CPU in stream: %fms\n", stepTime);
//		totalCommunicationTime += stepTime;

	}

	for (int i = 0; i < streamNum; i++)
	{
		cudaFree(dataBuf_d[i]);
		cudaFree(codeBuf_d[i]);
	}
	cudaFree(encodingMatrix_d);

	// record event and synchronize
	cudaEventRecord(totalStop);
	cudaEventSynchronize(totalStop);
	// get event elapsed time
	cudaEventElapsedTime(&totalTime, totalStart, totalStop);
//	printf("Total computation time: %fms\n", totalComputationTime);
//	printf("Total communication time: %fms\n", totalCommunicationTime);
	printf("Total GPU encoding time: %fms\n", totalTime);

	for (int i = 0; i < streamNum; i++)
	{
		cudaStreamDestroy(stream[i]);
	}

	if (id == 0)
	{
		char metadata_file_name[strlen(fileName) + 15];
		sprintf(metadata_file_name, "%s.METADATA", fileName);
		write_metadata(metadata_file_name, totalSize, parityBlockNum, nativeBlockNum, encodingMatrix);
	}
	// Pageable Host Memory is preferred here since the encodingMatrix is small
	free(encodingMatrix);
//	cudaFreeHost(encodingMatrix);
}

static void* GPU_thread_func(void * args)
{
	ThreadDataType* thread_data = (ThreadDataType *) args;
	cudaSetDevice(thread_data->id);
	struct timespec start, end;
	pthread_barrier_wait(&barrier);
	clock_gettime(CLOCK_REALTIME, &start);
	pthread_barrier_wait(&barrier);
	encode(thread_data->fileName, 
			thread_data->dataBuf, 
			thread_data->codeBuf, 
			thread_data->id, 
			thread_data->nativeBlockNum, 
			thread_data->parityBlockNum, 
			thread_data->chunkSize, 
			thread_data->totalSize, 
			thread_data->gridDimXSize, 
			thread_data->streamNum);
	pthread_barrier_wait(&barrier);
	clock_gettime(CLOCK_REALTIME, &end);
	if (thread_data->id == 0)
	{
		double totalTime = (double) (end.tv_sec - start.tv_sec) * 1000
				+ (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000L;
		printf("Total GPU encoding time using multiple devices: %fms\n", totalTime);
	}
	return NULL;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  encode_file
 *  Description:  encode the input file <fileName> with the given settings
 * =====================================================================================
 */
extern "C"
void encode_file(char *fileName, int nativeBlockNum, int parityBlockNum, int gridDimXSize, int streamNum)
{
	int chunkSize = 1;
	int totalSize;

	FILE *fp_in;
	FILE *fp_out;
	if ((fp_in = fopen(fileName, "rb")) == NULL)
	{
		printf("Cannot open source file!\n");
		exit(0);
	}

	fseek(fp_in, 0L, SEEK_END);
	// use ftell() to get the total size of the file
	totalSize = ftell(fp_in);
	chunkSize = (totalSize / nativeBlockNum) + (totalSize % nativeBlockNum != 0); 
//	chunkSize = (ftell(fp_in) / nativeBlockNum) + ( ftell(fp_in)%nativeBlockNum != 0 ); 
//	chunkSize = (int) (ceil( (long double) (ftell(fp_in) / nativeBlockNum)) ); 

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host
	int dataSize = nativeBlockNum * chunkSize * sizeof(uint8_t);
	int codeSize = parityBlockNum * chunkSize * sizeof(uint8_t);
	// Pinned host memory is expensive for allocation,
	// so pageable host memory is used here.
	dataBuf = (uint8_t*) malloc(dataSize);
	codeBuf = (uint8_t*) malloc(codeSize);
//	cudaMallocHost((void **)&dataBuf, dataSize);
//	memset(dataBuf, 0, dataSize);
//	cudaMallocHost((void **)&codeBuf, codeSize);
//	memset(codeBuf, 0, codeSize);
	
	for (int i = 0; i < nativeBlockNum; i++)
	{
		if (fseek(fp_in, i * chunkSize, SEEK_SET) == -1)
		{
			printf("fseek error!\n");
			exit(0);
		}

		if (fread(dataBuf + i * chunkSize, sizeof(uint8_t), chunkSize, fp_in) == EOF)
		{
			printf("fread error!\n");
			exit(0);
		}
	}
	fclose(fp_in);
	
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
	
	int GPU_num;
	cudaGetDeviceCount(&GPU_num);
	void* threads = malloc(GPU_num * sizeof(pthread_t));
	ThreadDataType* thread_data = (ThreadDataType *) malloc(GPU_num * sizeof(ThreadDataType));
	uint8_t *dataBufPerDevice[GPU_num];
	uint8_t *codeBufPerDevice[GPU_num];
	pthread_barrier_init(&barrier, NULL, GPU_num);
	int maxChunkSizePerDevice = (chunkSize / GPU_num) + (chunkSize % GPU_num != 0);
//	struct timespec start, end;
//	clock_gettime(CLOCK_REALTIME, &start);
	for (int i = 0; i < GPU_num; ++i)
	{
		thread_data[i].id = i;
		thread_data[i].nativeBlockNum = nativeBlockNum;
		thread_data[i].parityBlockNum = parityBlockNum;
		int deviceChunkSize = min(chunkSize - i * maxChunkSizePerDevice, maxChunkSizePerDevice);
		thread_data[i].chunkSize = deviceChunkSize;
		thread_data[i].totalSize = totalSize;
		thread_data[i].gridDimXSize = gridDimXSize;
		thread_data[i].streamNum = streamNum;
		thread_data[i].fileName = fileName;
		int deviceDataSize = nativeBlockNum * deviceChunkSize * sizeof(uint8_t);
		int deviceCodeSize = parityBlockNum * deviceChunkSize * sizeof(uint8_t);
		// Pageable Host Memory
//		dataBufPerDevice[i] = (uint8_t*) malloc(deviceDataSize);
//		codeBufPerDevice[i] = (uint8_t*) malloc(deviceCodeSize);
		// Pinned Host Memory
		cudaMallocHost((void **)&dataBufPerDevice[i], deviceDataSize);
		cudaMallocHost((void **)&codeBufPerDevice[i], deviceCodeSize);
		for (int j = 0; j < nativeBlockNum; ++j)
		{
			// Pageable Host Memory
//			memcpy(dataBufPerDevice[i] + j * deviceChunkSize, 
//							dataBuf + j * chunkSize + i * deviceChunkSize,
//							deviceChunkSize);
			// Pinned Host Memory
			cudaMemcpy(dataBufPerDevice[i] + j * deviceChunkSize, 
							dataBuf + j * chunkSize + i * deviceChunkSize,
							deviceChunkSize,
							cudaMemcpyHostToHost);
		}
		thread_data[i].dataBuf = dataBufPerDevice[i];
		thread_data[i].codeBuf = codeBufPerDevice[i];
		pthread_create(&((pthread_t*) threads)[i], NULL, GPU_thread_func, (void *) &thread_data[i]);
	}
	for (int i = 0; i < GPU_num; ++i)
	{
		pthread_join(((pthread_t*) threads)[i], NULL);
	}
//	clock_gettime(CLOCK_REALTIME, &end);
//	double totalTime = (double) (end.tv_sec - start.tv_sec) * 1000
//			+ (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000L;
//	printf("Total GPU encoding time using multiple devices: %fms\n", totalTime);
	for (int i = 0; i < GPU_num; ++i)
	{
		int deviceChunkSize = min(chunkSize - i * maxChunkSizePerDevice, maxChunkSizePerDevice);
		for (int j = 0; j < parityBlockNum; ++j)
		{
			// Pageable Host Memory
//			memcpy(codeBuf + j * chunkSize + i * deviceChunkSize,
//							codeBufPerDevice[i] + j * deviceChunkSize,
//							deviceChunkSize);
			// Pinned Host Memory
			cudaMemcpy(codeBuf + j * chunkSize + i * deviceChunkSize,
							codeBufPerDevice[i] + j * deviceChunkSize,
							deviceChunkSize,
							cudaMemcpyHostToHost);
		}
		// Pageable Host Memory
//		free(dataBufPerDevice[i]);
//		free(codeBufPerDevice[i]);
		// Pinned Host Memory
		cudaFreeHost(dataBufPerDevice[i]);
		cudaFreeHost(codeBufPerDevice[i]);
	}
	pthread_barrier_destroy(&barrier);
	cudaDeviceReset();

	char output_file_name[strlen(fileName) + 5];
	for (int i = 0; i < nativeBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_%s", i, fileName);
		if ((fp_out = fopen(output_file_name, "wb")) == NULL)
		{
			printf("Cannot open output file!\n");
			exit(0);
		}
		if (fwrite(dataBuf + i * chunkSize, sizeof(uint8_t), chunkSize, fp_out) != sizeof(uint8_t) * chunkSize)
		{
			printf("fwrite error!\n");
			exit(0);
		}
		fclose(fp_out);
	}
	for (int i = 0; i < parityBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_%s", i + nativeBlockNum, fileName);
		if ((fp_out = fopen(output_file_name, "wb")) == NULL)
		{
			printf("Cannot open output file!\n");
			exit(0);
		}
		if (fwrite(codeBuf + i * chunkSize, sizeof(uint8_t), chunkSize, fp_out) != sizeof(uint8_t)*chunkSize)
		{
			printf("fwrite error!\n");
			exit(0);
		}
		fclose(fp_out);
	}

	// Pinned host memory is expensive for deallocation,
	// so pageable host memory is used here.
	free(dataBuf);
	free(codeBuf);
//	cudaFreeHost(dataBuf);
//	cudaFreeHost(codeBuf);
}
