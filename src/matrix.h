/*
 * =====================================================================================
 *
 *       Filename:  matrix.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/21/2012 07:45:23 PM
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Shuai YUAN (yszheda AT gmail.com), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _MATRIX_H_
#define _MATRIX_H_

#define SQUARE_BLOCK_SIZE 16
#define SINGLE_BLOCK_SIZE 512

#define index(i, j, size) (((i) * (size)) + (j))

//#define TILE_WIDTH 2
#define TILE_WIDTH_ROW 2
#define TILE_WIDTH_COL 64
#define TILE_DEPTH 2

#define SINGLE_GRID_SIZE 16384 // MAX 
// #define SINGLE_GRID_SIZE 2147483647 // Max GridDim.x in K20 

// #define W 4
// #define NW (1 << W) /* In other words, NW equals 2 to the w-th power */

#define gf uint8_t
/*
__host__ __device__ uint8_t gf_add(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_sub(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp);
__host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp);
__host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power);
__host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power, uint8_t *gflog, uint8_t *gfexp);
// __device__ void matrix_mul(uint8_t *A, uint8_t *B, uint8_t *C, int n, int p, int m);
__global__ void matrix_mul(unsigned char *A, unsigned char *B, unsigned char *C, int n, int p, int m, int tileWidthRow, int tileWidthCol, int tileDepth);

// __global__ void encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize);
// __global__ void decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize);
*/
__global__ void gen_encoding_matrix(uint8_t *encodingMatrix, int row, int col);

__host__ float encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize, int tileWidthRow, int tileWidthCol, int tileDepth);

__host__ float decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize, int tileWidthRow, int tileWidthCol, int tileDepth);

extern "C"
void invert_matrix(uint8_t *matrix_dev, uint8_t *result_dev, int size);
#endif
