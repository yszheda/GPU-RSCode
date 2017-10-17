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

#include "encode.h"
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
    uint8_t* encodingMatrix;
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
void encode(uint8_t *dataBuf, uint8_t *codeBuf, uint8_t *encodingMatrix, int id, int nativeBlockNum, int parityBlockNum, int chunkSize, int totalSize, int gridDimXSize, int streamNum)
{
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

    uint8_t *encodingMatrix_d;	//device
    int matrixSize = parityBlockNum * nativeBlockNum * sizeof(uint8_t);
    cudaMalloc((void **)&encodingMatrix_d, matrixSize);

    // record event
    cudaEventRecord(stepStart);
    const int maxBlockDimSize = 16;
    int blockDimX = min(parityBlockNum, maxBlockDimSize);
    int blockDimY = min(nativeBlockNum, maxBlockDimSize);
    int gridDimX = (int) ceil((float) parityBlockNum / blockDimX);
    int gridDimY = (int) ceil((float) nativeBlockNum / blockDimY);
    dim3 grid(gridDimX, gridDimY);
    dim3 block(blockDimX, blockDimY);
    gen_encoding_matrix<<<grid, block>>>(encodingMatrix_d, parityBlockNum, nativeBlockNum);
    // cudaDeviceSynchronize();
    // record event and synchronize
    cudaEventRecord(stepStop);
    cudaEventSynchronize(stepStop);
    // get event elapsed time
    cudaEventElapsedTime(&stepTime, stepStart, stepStop);
    printf("Device%d: Generating encoding matrix completed: %fms\n", id, stepTime);
    totalComputationTime += stepTime;

    if (id == 0)
    {
        // record event
        cudaEventRecord(stepStart);
        cudaMemcpy(encodingMatrix, encodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
        // record event and synchronize
        cudaEventRecord(stepStop);
        cudaEventSynchronize(stepStop);
        // get event elapsed time
        cudaEventElapsedTime(&stepTime, stepStart, stepStop);
        printf("Device%d: Copy encoding matrix from GPU to CPU: %fms\n", id, stepTime);
        totalCommunicationTime += stepTime;
    }

    // Use cuda stream to encode the file
    // to achieve computation and comunication overlapping
    // Use DFS way
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
        for (int j = 0; j < nativeBlockNum; j++)
        {
            cudaMemcpyAsync(dataBuf_d[i] + j * streamChunkSize,
                    dataBuf + j * chunkSize + i * streamMinChunkSize,
                    streamChunkSize * sizeof(uint8_t),
                    cudaMemcpyHostToDevice,
                    stream[i]);
        }

        stepTime = encode_chunk(dataBuf_d[i], encodingMatrix_d, codeBuf_d[i], nativeBlockNum, parityBlockNum, streamChunkSize, gridDimXSize, stream[i]);

        for (int j = 0; j < parityBlockNum; j++)
        {
            cudaMemcpyAsync(codeBuf + j * chunkSize + i * streamMinChunkSize,
                    codeBuf_d[i] + j * streamChunkSize,
                    streamChunkSize * sizeof(uint8_t),
                    cudaMemcpyDeviceToHost,
                    stream[i]);
        }
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
    printf("Device%d: Total GPU encoding time: %fms\n", id, totalTime);

    for (int i = 0; i < streamNum; i++)
    {
        cudaStreamDestroy(stream[i]);
    }
}

static void* GPU_thread_func(void * args)
{
    ThreadDataType* thread_data = (ThreadDataType *) args;
    cudaSetDevice(thread_data->id);

    uint8_t *encodingMatrix;
    int parityBlockNum = thread_data->parityBlockNum;
    int nativeBlockNum = thread_data->nativeBlockNum;
    int matrixSize = parityBlockNum * nativeBlockNum * sizeof(uint8_t);
    // NOTE: Pageable Host Memory is preferred here since the encodingMatrix is small
    // cudaMallocHost((void **)&encodingMatrix, matrixSize);
    encodingMatrix = (uint8_t*) malloc(matrixSize);
    thread_data->encodingMatrix = encodingMatrix;

    struct timespec start, end;
    pthread_barrier_wait(&barrier);
    clock_gettime(CLOCK_REALTIME, &start);
    pthread_barrier_wait(&barrier);
    encode(thread_data->dataBuf,
            thread_data->codeBuf,
            thread_data->encodingMatrix,
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

    if (thread_data->id == 0)
    {
        char *fileName = thread_data->fileName;
        int totalSize = thread_data->totalSize;
        char metadata_file_name[strlen(fileName) + 15];
        sprintf(metadata_file_name, "%s.METADATA", fileName);
        write_metadata(metadata_file_name, totalSize, parityBlockNum, nativeBlockNum, encodingMatrix);
    }
    // NOTE: Pageable Host Memory is preferred here since the encodingMatrix is small
    free(encodingMatrix);
    // cudaFreeHost(encodingMatrix);

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
    chunkSize = (totalSize + nativeBlockNum - 1) / nativeBlockNum;

    uint8_t *dataBuf;		//host
    uint8_t *codeBuf;		//host
    int dataSize = nativeBlockNum * chunkSize * sizeof(uint8_t);
    int codeSize = parityBlockNum * chunkSize * sizeof(uint8_t);
    // NOTE: Pinned host memory is expensive for allocation,
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
    int maxGridDimXSize = min(deviceProperties.maxGridSize[0], deviceProperties.maxGridSize[1]);
    if (gridDimXSize > maxGridDimXSize || gridDimXSize <= 0)
    {
        printf("Valid grid size: (0, %d]\n", maxGridDimXSize);
        gridDimXSize = maxGridDimXSize;
    }

    int GPU_num;
    cudaGetDeviceCount(&GPU_num);

    void* threads = malloc(GPU_num * sizeof(pthread_t));
    ThreadDataType* thread_data = (ThreadDataType *) malloc(GPU_num * sizeof(ThreadDataType));

    uint8_t *dataBufPerDevice[GPU_num];
    uint8_t *codeBufPerDevice[GPU_num];

    pthread_barrier_init(&barrier, NULL, GPU_num);

    int minChunkSizePerDevice = chunkSize / GPU_num;
    for (int i = 0; i < GPU_num; ++i)
    {
        cudaSetDevice(i);

        thread_data[i].id = i;
        thread_data[i].nativeBlockNum = nativeBlockNum;
        thread_data[i].parityBlockNum = parityBlockNum;
        int deviceChunkSize = minChunkSizePerDevice;
        if (i == GPU_num - 1)
        {
            deviceChunkSize = chunkSize - i * minChunkSizePerDevice;
        }
        thread_data[i].chunkSize = deviceChunkSize;
        thread_data[i].totalSize = totalSize;
        thread_data[i].gridDimXSize = gridDimXSize;
        thread_data[i].streamNum = streamNum;
        thread_data[i].fileName = fileName;

        int deviceDataSize = nativeBlockNum * deviceChunkSize * sizeof(uint8_t);
        int deviceCodeSize = parityBlockNum * deviceChunkSize * sizeof(uint8_t);
        cudaMallocHost((void **)&dataBufPerDevice[i], deviceDataSize);
        cudaMallocHost((void **)&codeBufPerDevice[i], deviceCodeSize);
        for (int j = 0; j < nativeBlockNum; ++j)
        {
            // Pinned Host Memory
            cudaMemcpy(dataBufPerDevice[i] + j * deviceChunkSize,
                    dataBuf + j * chunkSize + i * minChunkSizePerDevice,
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

    for (int i = 0; i < GPU_num; ++i)
    {
        int deviceChunkSize = minChunkSizePerDevice;
        if (i == GPU_num - 1) {
            deviceChunkSize = chunkSize - i * minChunkSizePerDevice;
        }

        for (int j = 0; j < parityBlockNum; ++j)
        {
            // Pinned Host Memory
            cudaMemcpy(codeBuf + j * chunkSize + i * minChunkSizePerDevice,
                    codeBufPerDevice[i] + j * deviceChunkSize,
                    deviceChunkSize,
                    cudaMemcpyHostToHost);
        }

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

    // NOTE: Pinned host memory is expensive for deallocation,
    // so pageable host memory is used here.
    free(dataBuf);
    free(codeBuf);
    // cudaFreeHost(dataBuf);
    // cudaFreeHost(codeBuf);
}
