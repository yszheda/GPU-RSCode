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
#ifndef _DECODE_H_
#define _DECODE_H_


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  decode_file
 *  Description:  decode the original input file <inFile> with the given settings
 *  @param inFile: name of the origin input file
 *  @param confFile: name of the configuration file
 *  @param outFile: name of the output file
 *  @param gridDimXSize: maximum grid size of X dimension
 *  @param streamNum: number of CUDA streams
 *
 *  The configuration file contains which files are used for decoding.
 * =====================================================================================
 */
#ifdef __cplusplus
extern "C" {
#endif
void decode_file(char *inFile, char *confFile, char *outFile, int gridDimXSize, int streamNum);
#ifdef __cplusplus
}
#endif


#endif
