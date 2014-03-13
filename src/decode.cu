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
#include <pthread.h>
#include "matrix.h"
#include "cpu-decode.h"
extern "C"		
void CPU_invert_matrix(uint8_t *matrix, uint8_t *result, int size);


// #define DEBUG

struct ThreadDataType {
	int id;
	int nativeBlockNum;
	int parityBlockNum;
	int chunkSize;
	int totalSize;
	uint8_t* dataBuf;
	uint8_t* codeBuf;
	uint8_t* decodingMatrix;
};	/* ----------  end of struct ThreadDataType  ---------- */

typedef struct ThreadDataType ThreadDataType;

static pthread_barrier_t barrier;

void show_squre_matrix(uint8_t *matrix, int size)
{
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			printf("%d ", matrix[i*size+j]);
		}
		printf("\n");
	}
}

void copy_matrix(uint8_t *src, uint8_t *des, int srcRowIndex, int desRowIndex, int rowSize)
{
	for(int i = 0; i < rowSize; i++)
	{
		des[desRowIndex * rowSize + i] = src[srcRowIndex * rowSize + i];
	}
}

extern "C"
void decode(uint8_t *dataBuf, uint8_t *codeBuf, uint8_t *decodingMatrix, int id, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
//	cudaSetDevice(id);

	int dataSize = nativeBlockNum * chunkSize * sizeof(uint8_t);
	int codeSize = nativeBlockNum * chunkSize * sizeof(uint8_t);
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

	cudaMalloc((void **) &dataBuf_d, dataSize);
//	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc((void **) &codeBuf_d, codeSize);
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

	int matrixSize = nativeBlockNum * nativeBlockNum * sizeof(uint8_t);
//	uint8_t *encodingMatrix_d;	//device
	uint8_t *decodingMatrix_d;	//device
//	cudaMalloc((void **) &encodingMatrix_d, matrixSize);
	cudaMalloc((void **) &decodingMatrix_d, matrixSize);

	// record event
	cudaEventRecord(stepStart);
	cudaMemcpy(decodingMatrix_d, decodingMatrix, matrixSize, cudaMemcpyHostToDevice);
	// record event and synchronize
	cudaEventRecord(stepStop);
	cudaEventSynchronize(stepStop);
	// get event elapsed time
	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
	printf("Copy decoding matrix from CPU to GPU: %fms\n", stepTime);
	totalCommunicationTime += stepTime;

//	// record event
//	cudaEventRecord(stepStart);
//	invert_matrix(encodingMatrix_d, decodingMatrix_d, nativeBlockNum);
//	// record event and synchronize
//	cudaEventRecord(stepStop);
//	cudaEventSynchronize(stepStop);
//	// get event elapsed time
//	cudaEventElapsedTime(&stepTime, stepStart, stepStop);
//	printf("Generating decoding matrix completed: %fms\n", stepTime);
//	totalComputationTime += stepTime;
//
//#ifdef DEBUG
//	uint8_t *decodingMatrix;	//host
//	decodingMatrix = (uint8_t*) malloc(matrixSize);
//	cudaMemcpy(decodingMatrix, decodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
//	show_squre_matrix(decodingMatrix, nativeBlockNum);
//	free(decodingMatrix);
//#endif
//
	stepTime = decode_chunk(dataBuf_d, decodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
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

static void* GPU_thread_func(void * args)
{
	ThreadDataType* thread_data = (ThreadDataType *) args;
	cudaSetDevice(thread_data->id);
	struct timespec start, end;
	pthread_barrier_wait(&barrier);
	clock_gettime(CLOCK_REALTIME, &start);
	pthread_barrier_wait(&barrier);
	decode(thread_data->dataBuf, 
			thread_data->codeBuf, 
			thread_data->decodingMatrix, 
			thread_data->id, 
			thread_data->nativeBlockNum, 
			thread_data->parityBlockNum, 
			thread_data->chunkSize);
	pthread_barrier_wait(&barrier);
	clock_gettime(CLOCK_REALTIME, &end);
	if (thread_data->id == 0)
	{
		double totalTime = (double) (end.tv_sec - start.tv_sec) * 1000
				+ (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000L;
		printf("Total GPU decoding time using multiple devices: %fms\n", totalTime);
	}
	return NULL;
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
		printf("Cannot open metadata file!\n");
		exit(0);
	}
	fscanf(fp_meta, "%d", &totalSize);
	fscanf(fp_meta, "%d %d", &parityBlockNum, &nativeBlockNum);

	chunkSize = (totalSize / nativeBlockNum) + (totalSize%nativeBlockNum != 0); 
//	chunkSize = (int) (ceil((double) totalSize / nativeBlockNum)); 
#ifdef DEBUG
printf("chunk size: %d\n", chunkSize);
#endif

	totalMatrixSize = nativeBlockNum * (nativeBlockNum + parityBlockNum);
	totalEncodingMatrix = (uint8_t*) malloc(totalMatrixSize);
	matrixSize = nativeBlockNum * nativeBlockNum;
	encodingMatrix = (uint8_t*) malloc(matrixSize);
	for (int i = 0; i < totalMatrixSize; ++i)
	{
//		fscanf(fp_meta, "%d", &totalEncodingMatrix[i]);
		int j;
		fscanf(fp_meta, "%d", &j);
		totalEncodingMatrix[i] = (uint8_t) j; 
	}
	fclose(fp_meta);

	dataSize = nativeBlockNum * chunkSize * sizeof(uint8_t);
	codeSize = nativeBlockNum * chunkSize * sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc(dataSize);
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc(codeSize);
	memset(codeBuf, 0, codeSize);

	FILE *fp_conf;
	char input_file_name[strlen(inFile) + 20];
	int index;
	if((fp_conf = fopen(confFile, "r")) == NULL)
	{
		printf("Cannot open configuration file!\n");
		exit(0);
	}

	for(int i = 0; i < nativeBlockNum; i++)
	{
		fscanf(fp_conf, "%s", input_file_name);
		index = atoi(input_file_name + 1);

		copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

		if((fp_in = fopen(input_file_name, "rb")) == NULL)
		{
			printf("Cannot open input file %s!\n", input_file_name);
			exit(0);
		}
		fseek(fp_in, 0L, SEEK_SET);
		// this part can be process in parallel with computing inversed matrix
		fread(codeBuf + i * chunkSize, sizeof(uint8_t), chunkSize, fp_in);
		fclose(fp_in);
	}
	fclose(fp_conf);
	
	uint8_t *decodingMatrix;
	decodingMatrix = (uint8_t*) malloc(matrixSize);
    CPU_invert_matrix(encodingMatrix, decodingMatrix, nativeBlockNum);

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
		int deviceDataSize = nativeBlockNum * deviceChunkSize * sizeof(uint8_t);
		int deviceCodeSize = nativeBlockNum * deviceChunkSize * sizeof(uint8_t);
//		dataBufPerDevice[i] = (uint8_t*) malloc(deviceDataSize);
		cudaMallocHost((void **)&dataBufPerDevice[i], deviceDataSize);
//		codeBufPerDevice[i] = (uint8_t*) malloc(deviceCodeSize);
		cudaMallocHost((void **)&codeBufPerDevice[i], deviceCodeSize);
		for (int j = 0; j < nativeBlockNum; ++j)
		{
//			memcpy(codeBufPerDevice[i] + j * deviceChunkSize, 
//							codeBuf + j * chunkSize + i * deviceChunkSize,
//							deviceChunkSize);
			cudaMemcpy(codeBufPerDevice[i] + j * deviceChunkSize, 
							codeBuf + j * chunkSize + i * deviceChunkSize,
							deviceChunkSize,
							cudaMemcpyHostToHost);
		}
		thread_data[i].dataBuf = dataBufPerDevice[i];
		thread_data[i].codeBuf = codeBufPerDevice[i];
		thread_data[i].decodingMatrix = decodingMatrix;
		pthread_create(&((pthread_t*) threads)[i], NULL, GPU_thread_func, (void *) &thread_data[i]);
	}
	for (int i = 0; i < GPU_num; ++i)
	{
		pthread_join(((pthread_t*) threads)[i], NULL);
	}
//	clock_gettime(CLOCK_REALTIME, &end);
//	double totalTime = (double) (end.tv_sec - start.tv_sec) * 1000
//			+ (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000L;
//	printf("Total GPU decoding time using multiple devices: %fms\n", totalTime);
	for (int i = 0; i < GPU_num; ++i)
	{
		int deviceChunkSize = min(chunkSize - i * maxChunkSizePerDevice, maxChunkSizePerDevice);
		for (int j = 0; j < nativeBlockNum; ++j)
		{
//			memcpy(dataBuf + j * chunkSize + i * deviceChunkSize,
//							dataBufPerDevice[i] + j * deviceChunkSize,
//							deviceChunkSize);
			cudaMemcpy(dataBuf + j * chunkSize + i * deviceChunkSize,
							dataBufPerDevice[i] + j * deviceChunkSize,
							deviceChunkSize,
							cudaMemcpyHostToHost);
		}
//		free(dataBufPerDevice[i]);
		cudaFreeHost(dataBufPerDevice[i]);
//		free(codeBufPerDevice[i]);
		cudaFreeHost(codeBufPerDevice[i]);
	}
	pthread_barrier_destroy(&barrier);
	cudaDeviceReset();

	if(outFile == NULL)
	{
		if((fp_out = fopen(inFile, "wb")) == NULL)
		{
			printf("Cannot open output file %s!\n", inFile);
			exit(0);
		}
	}
	else
	{
		if((fp_out = fopen(outFile, "wb")) == NULL)
		{
			printf("Cannot open output file %s!\n", outFile);
			exit(0);
		}
	}
	fwrite(dataBuf, sizeof(uint8_t), totalSize, fp_out);
	fclose(fp_out);

	free(dataBuf);
	free(codeBuf);
}
