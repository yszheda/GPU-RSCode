/*
 * =====================================================================================
 *
 *       Filename:  Matrix.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/05/13 17:24:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuai YUAN(yszheda AT gmail.com), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef  MATRIX_H
#define  MATRIX_H
#include <iostream>

template <typename T>
void matrixMul ( const int rows_lhs, const int cols_lhs, const int cols_rhs, const T* lhs, const T* rhs, T* result );

template <typename T>
void matrixInvert ( const int size, T* matrix, T* result );

#endif   /* ----- #ifndef MATRIX_H  ----- */
