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

#include <cuda.h>
#include <stdint.h>

#define SQUARE_BLOCK_SIZE 16
#define SINGLE_BLOCK_SIZE 512

#define index(i, j, size) (((i) * (size)) + (j))

//#define TILE_WIDTH 2
#define TILE_WIDTH_ROW 2
#define TILE_WIDTH_COL 64
#define TILE_DEPTH 2

#define SINGLE_GRID_SIZE 16384 // MAX
// #define SINGLE_GRID_SIZE 2147483647 // Max GridDim.x in K20

#define gf uint8_t

// AlignType: the word-alignment length
// typedef unsigned int AlignType;
typedef long long AlignType;

// __host__ __device__ uint8_t gf_add(uint8_t a, uint8_t b);
// __host__ __device__ uint8_t gf_sub(uint8_t a, uint8_t b);
// __host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b);
// __host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp);
// __host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b);
// __host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp);
// __host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power);
// __host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power, uint8_t *gflog, uint8_t *gfexp);
// __global__ void matrix_mul(unsigned char *A, unsigned char *B, unsigned char *C, int n, int p, int m, int tileWidthRow, int tileWidthCol, int tileDepth);

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  gen_encoding_matrix
 *  Description:  generate encoding matrix
 *  @param encodingMatrix: the pointer to the encoding matrix
 *  @param parityBlockNum: number of parity code chunks
 *  @param nativeBlockNum: number of native data chunks
 * =====================================================================================
 */
__global__ void gen_encoding_matrix(uint8_t *encodingMatrix, int parityBlockNum, int nativeBlockNum);

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  encode_chunk
 *  Description:  encode the given buffer of data chunks
 *  @param dataChunk: pointer to the buffer of data chunks
 *  @param parityCoeff: parity coefficient matrix (encoding matrix)
 *  @param codeChunk: pointer to the buffer of code chunks
 *  @param nativeBlockNum: number of native data chunks
 *  @param parityBlockNum: number of parity code chunks
 *  @param chunkSize: size of each chunk
 *  @param gridDimXSize: maximum grid size of X dimension
 *  @param streamID: ID of the current CUDA stream
 *  @return the kernel execution time
 * =====================================================================================
 */
__host__ float encode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize, int gridDimXSize, cudaStream_t streamID);

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  decode_chunk
 *  Description:  decode the given buffer of code chunks
 *  @param dataChunk: pointer to the buffer of data chunks
 *  @param parityCoeff: parity coefficient matrix (encoding matrix)
 *  @param codeChunk: pointer to the buffer of code chunks
 *  @param nativeBlockNum: number of native data chunks
 *  @param parityBlockNum: number of parity code chunks
 *  @param chunkSize: size of each chunk
 *  @param gridDimXSize: maximum grid size of X dimension
 *  @param streamID: ID of the current CUDA stream
 *  @return the kernel execution time
 * =====================================================================================
 */
__host__ float decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize, int gridDimXSize, cudaStream_t streamID);

#ifdef __cplusplus
extern "C" {
#endif
void GPU_invert_matrix(uint8_t *matrix_dev, uint8_t *result_dev, int size);
#ifdef __cplusplus
}
#endif

#endif
