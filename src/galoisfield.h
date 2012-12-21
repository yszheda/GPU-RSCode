/*
 * =====================================================================================
 *
 *       Filename:  galoisfield.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/27/2012 10:46:18 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _GALOISFIELD_H_
#define _GALOISFIELD_H_

__shared__ uint8_t *gflog;
__shared__ uint8_t *gfexp;

__host__ __device__ int setup_tables(int w);
__host__ __device__ uint8_t gf_add(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_sub(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_mul(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp);
__host__ __device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp);
__host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_div(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp);
__host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t b);
__host__ __device__ uint8_t gf_pow(uint8_t a, uint8_t power, uint8_t *gflog, uint8_t *gfexp);

#endif
