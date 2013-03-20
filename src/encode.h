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
#ifdef _ENCODE_H_
#define _ENCODE_H_

// extern "C"
// void encode(uint8_t *dataBuf, uint8_t *codeBuf, uint8_t *encodingMatrix, int nativeBlockNum, int parityBlockNum, int chunkSize, int totalSize);

extern "C"
void encode_file(char *file, int nativeBlockNum, int parityBlockNum);
#endif
