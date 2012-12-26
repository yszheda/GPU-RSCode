#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <stdint.h>
#include "matrix.h"

//#define SQUARE_BLOCK_SIZE 16    // MAX 
//#define LINEAR_BLOCK_SIZE 512   // MAX 
//
//#define DISPLAY_SETTINGS false
//#define DISPLAY false
//
////#define IDC2D(i,j,ld) (((j)*(ld))+(i))
//#define IDC2D(i,j,ld) (((i)*(ld))+(j))
//
//#define W 8
//#define NW (1 << W) /* In other words, NW equals 2 to the w-th power */

//#define TILE_WIDTH 4
//#define TILE_DEPTH 4

//#define BUFFER_SIZE 256

#define DEBUG

//__global__ void decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
//{
//	matrix_mul(parityCoeff, codeChunk, dataChunk, nativeBlockNum, nativeBlockNum, chunkSize);
//}

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
	cudaMalloc( (void **)&dataBuf_d, dataSize );
//	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc( (void **)&codeBuf_d, codeSize );
//	cudaMemset(codeBuf_d, 0, codeSize);
	cudaMemcpy(codeBuf_d, codeBuf, codeSize, cudaMemcpyHostToDevice);

	int matrixSize = nativeBlockNum * nativeBlockNum;
	uint8_t *encodingMatrix_d;	//device
	uint8_t *decodingMatrix_d;	//device
	cudaMalloc( (void **)&encodingMatrix_d, matrixSize );
	cudaMalloc( (void **)&decodingMatrix_d, matrixSize );
	cudaMemcpy(encodingMatrix_d, encodingMatrix, matrixSize, cudaMemcpyHostToDevice);

    invert_matrix(encodingMatrix_d, decodingMatrix_d, nativeBlockNum);

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
	int gridDimX = (int)( ceil((float)chunkSize / TILE_WIDTH_COL) );
	int gridDimY = (int)( ceil((float)nativeBlockNum / TILE_WIDTH_ROW) );
	dim3 grid(gridDimX, gridDimY);
	dim3 block(TILE_WIDTH_ROW, TILE_WIDTH_COL);
	decode_chunk<<<grid, block>>>(dataBuf_d, decodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);

	cudaMemcpy(dataBuf, dataBuf_d, dataSize, cudaMemcpyDeviceToHost);

	cudaFree(decodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);
}

extern "C"
void decode_file(char *confFile, int nativeBlockNum, int parityBlockNum)
{
	int chunkSize = 1;
	int totalSize;

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host

	int dataSize;
	int codeSize;

	FILE *fp_in;
	FILE *fp_out;

	int totalMatrixSize;
	int matrixSize;
	uint8_t *totalEncodingMatrix;	//host
	uint8_t *encodingMatrix;	//host
	if( ( fp_in = fopen(".METADATA","rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}
	fscanf(fp_in, "%d", &totalSize);
	fscanf(fp_in, "%d %d", &parityBlockNum, &nativeBlockNum);
	chunkSize = (int) (ceil( (float)totalSize / nativeBlockNum )); 
#ifdef DEBUG
printf("chunk size: %d\n", chunkSize);
#endif
	totalMatrixSize = nativeBlockNum * ( nativeBlockNum + parityBlockNum );
	totalEncodingMatrix = (uint8_t*) malloc( totalMatrixSize );
	matrixSize = nativeBlockNum * nativeBlockNum;
	encodingMatrix = (uint8_t*) malloc( matrixSize );
	for(int i =0; i<nativeBlockNum*(nativeBlockNum+parityBlockNum); i++)
	{
		fscanf(fp_in, "%d", totalEncodingMatrix+i);
	}

	dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( dataSize );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( codeSize);
	memset(codeBuf, 0, codeSize);

	if(confFile != NULL)
	{
		FILE *fp_conf;
		char input_file_name[20];
		int index;
		fp_conf = fopen(confFile, "r");

		for(int i=0; i<nativeBlockNum; i++)
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
	}
	else
	{
		for(int i=0; i<nativeBlockNum; i++)
		{
			char input_file_name[20];
			int index;
			printf("Please enter the file name of fragment:\n");
			scanf("%s", input_file_name);
			index = atoi(input_file_name+1);
			printf("#%dth fragment\n", index);

			copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

			fp_in = fopen(input_file_name, "rb");
			fseek(fp_in, 0L, SEEK_SET);
			// this part can be process in parallel with computing inversed matrix
			fread(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
			fclose(fp_in);

		}
	}
/*
	for(int i=0; i<nativeBlockNum; i++)
	{
		char input_file_name[20];
		int index;
		printf("Please enter the file name of fragment:\n");
		scanf("%s", input_file_name);
		index = atoi(input_file_name+1);
		printf("#%dth fragment\n", index);

		copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

		fp_in = fopen(input_file_name, "rb");
		fseek(fp_in, 0L, SEEK_SET);
		// this part can be process in parallel with computing inversed matrix
		fread(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
		fclose(fp_in);

	}
*/
	
	decode(dataBuf, codeBuf, encodingMatrix, nativeBlockNum, parityBlockNum, chunkSize);

	char output_file_name[20];
	printf("Enter the name of the decoded file:\n");
	scanf("%s", output_file_name);
	fp_out = fopen(output_file_name, "wb");
	fwrite(dataBuf, sizeof(uint8_t), totalSize, fp_out);
	fclose(fp_out);

	free(dataBuf);
	free(codeBuf);

}

/*
int main()
{
	int nativeBlockNum = 4;
	int parityBlockNum = 2;
	int chunkSize = 1;
	int totalSize;

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host
	uint8_t *dataBuf_d;		//device
	uint8_t *codeBuf_d;		//device

	int dataSize;
	int codeSize;

	FILE *fp_in;
	FILE *fp_out;

	int totalMatrixSize;
	int matrixSize;
	uint8_t *totalEncodingMatrix;	//host
	uint8_t *encodingMatrix;	//host
	if( ( fp_in = fopen(".METADATA","rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}
	fscanf(fp_in, "%d", &totalSize);
	fscanf(fp_in, "%d %d", &parityBlockNum, &nativeBlockNum);
	chunkSize = (int) (ceil( (float)totalSize / nativeBlockNum )); 
	totalMatrixSize = nativeBlockNum * ( nativeBlockNum + parityBlockNum );
	totalEncodingMatrix = (uint8_t*) malloc( totalMatrixSize );
	matrixSize = nativeBlockNum * nativeBlockNum;
	encodingMatrix = (uint8_t*) malloc( matrixSize );
	for(int i =0; i<nativeBlockNum*(nativeBlockNum+parityBlockNum); i++)
	{
		fscanf(fp_in, "%d", totalEncodingMatrix+i);
	}

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

	for(int i=0; i<nativeBlockNum; i++)
	{
		char input_file_name[20];
		int index;
		printf("Please enter the file name of fragment:\n");
		scanf("%s", input_file_name);
		index = atoi(input_file_name+1);
		printf("#%dth fragment\n", index);

		copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

		fp_in = fopen(input_file_name, "rb");
		fseek(fp_in, 0L, SEEK_SET);
		// this part can be process in parallel with computing inversed matrix
		fread(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
		fclose(fp_in);

	}
	cudaMemcpy(codeBuf_d, codeBuf, codeSize, cudaMemcpyHostToDevice);

#ifdef DEBUG
	show_squre_matrix(encodingMatrix, nativeBlockNum);
#endif

//	if( ( fp_in = fopen("native_2","rb") ) == NULL )
//	{
//		printf("Can not open source file!\n");
//		exit(0);
//	}
//	fseek(fp_in, 0L, SEEK_END);
//	chunkSize = ftell(fp_in);
#ifdef DEBUG
//	int matrixSize = nativeBlockNum*nativeBlockNum*sizeof(uint8_t);
	uint8_t testMatrix[4][4] = {{0, 0, 1, 0}, {0, 1, 0, 0}, {1, 1, 1, 1}, {1, 2, 3, 4}};
//	uint8_t testMatrix[16] = {0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 3, 4};
//	uint8_t testMatrix[16] = {1, 1, 1, 1, 1, 2, 3, 4, 0, 0, 1, 0, 0, 0, 0, 1};
	cudaMemcpy(encodingMatrix_d, testMatrix, matrixSize, cudaMemcpyHostToDevice);
#endif

	uint8_t *encodingMatrix_d;	//device
	uint8_t *decodingMatrix;	//host
	uint8_t *decodingMatrix_d;	//device
	decodingMatrix = (uint8_t*) malloc( matrixSize );
	cudaMalloc( (void **)&encodingMatrix_d, matrixSize );
	cudaMalloc( (void **)&decodingMatrix_d, matrixSize );
	cudaMemcpy(encodingMatrix_d, encodingMatrix, matrixSize, cudaMemcpyHostToDevice);

    invert_matrix( encodingMatrix_d, decodingMatrix_d, nativeBlockNum );

	cudaMemcpy(decodingMatrix, decodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
	show_squre_matrix(decodingMatrix, nativeBlockNum);

//	fseek(fp_in, 0L, SEEK_SET);
//	fread(codeBuf, sizeof(uint8_t), chunkSize, fp_in);
//	fclose(fp_in);
//
//	fp_in = fopen("native_1", "rb");
//	fseek(fp_in, 0L, SEEK_SET);
//	fread(codeBuf+chunkSize, sizeof(uint8_t), chunkSize, fp_in);
//	fclose(fp_in);
//
//	fp_in = fopen("code_0", "rb");
//	fseek(fp_in, 0L, SEEK_SET);
//	fread(codeBuf+2*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
//	fclose(fp_in);
//
//	fp_in = fopen("code_1", "rb");
//	fseek(fp_in, 0L, SEEK_SET);
//	fread(codeBuf+3*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
//	fclose(fp_in);


	int gridDimX = (int)(ceil((float)chunkSize/TILE_WIDTH));
	int gridDimY = (int)(ceil((float)parityBlockNum/TILE_WIDTH));
	dim3 grid(gridDimX, gridDimY);
	dim3 block(TILE_WIDTH, TILE_WIDTH);
	decode_chunk<<<grid, block>>>(dataBuf_d, decodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
	cudaMemcpy(dataBuf, dataBuf_d, dataSize, cudaMemcpyDeviceToHost);

	char output_file_name[20];
	printf("Enter the name of the decoded file:\n");
	scanf("%s", output_file_name);
	fp_out = fopen(output_file_name, "wb");
	fwrite(dataBuf, sizeof(uint8_t), totalSize, fp_out);
	fclose(fp_out);

	cudaFree(decodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);

	free(decodingMatrix);
	free(dataBuf);
	free(codeBuf);

}
*/
