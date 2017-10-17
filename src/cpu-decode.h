/*
 * =====================================================================================
 *
 *       Filename:  cpu-decode.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/04/14 09:36:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuai YUAN(yszheda AT gmail.com), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _CPU_DECODE_H_
#define _CPU_DECODE_H_

#define index(i, j, size) (((i) * (size)) + (j))

#ifdef __cplusplus
extern "C" {
#endif
    void CPU_invert_matrix(uint8_t *matrix, uint8_t *result, int size);
#ifdef __cplusplus
}
#endif

#endif
