/*
 * =====================================================================================
 *
 *       Filename:  decode.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/25/2012 04:35:32 PM
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Shuai YUAN (yszheda AT gmail.com), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifdef _DECODE_H_
#define _DECODE_H_

extern "C"
void decode_file(char *file, int nativeBlockNum, int parityBlockNum, int tileWidthRow, int tileWidthCol, int tileDepth);

extern "C"
void decode(uint8_t *dataBuf, uint8_t *codeBuf, uint8_t *encodingMatrix, int nativeBlockNum, int parityBlockNum, int chunkSize, int tileWidthRow, int tileWidthCol, int tileDepth);


#endif
