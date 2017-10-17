/*
 * =====================================================================================
 *
 *       Filename:  encode.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/25/2012 04:33:06 PM
 *       Revision:  none
 *       Compiler:  nvcc 
 *
 *         Author:  Shuai YUAN (yszheda AT gmail.com), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef _ENCODE_H_
#define _ENCODE_H_

// extern "C"
// void encode(char *fileName, uint8_t *dataBuf, uint8_t *codeBuf, int nativeBlockNum, int parityBlockNum, int chunkSize, int totalSize, int gridDimXSize, int streamNum);

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  encode_file
 *  Description:  encode the input file <fileName> with the given settings
 *  @param fileName: name of the input file
 *  @param nativeBlockNum: number of native data chunks
 *  @param parityBlockNum: number of parity code chunks
 *  @param gridDimXSize: maximum grid size of X dimension
 *  @param streamNum: number of CUDA streams
 * =====================================================================================
 */
#ifdef __cplusplus
extern "C" {
#endif
void encode_file(char *fileName, int nativeBlockNum, int parityBlockNum, int gridDimXSize, int streamNum);
#ifdef __cplusplus
}
#endif

#endif
