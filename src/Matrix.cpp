/*
 * =====================================================================================
 *
 *       Filename:  Matrix.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/05/13 16:22:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuai YUAN(yszheda AT gmail.com), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <iostream>
#include "Matrix.hpp"

inline int index(int row, int col, int table_size)
{
	return row * table_size + col;
}

template <typename T>
void matrixMul ( const int rows_lhs, const int cols_lhs, const int cols_rhs, const T* lhs, const T* rhs, T* result )
{
	for (int i = 0; i < rows_lhs; ++i)
	{
		for (int j = 0; j < cols_rhs; ++j)
		{
			for (int k = 0; k < cols_lhs; ++k)
			{
				result[ index(i, j, cols_rhs) ] += lhs[ index(i, k, cols_lhs) ] * rhs[ index(k, j, cols_rhs) ];
			}
		}
	}
}		/* -----  end of template function matrixMul  ----- */


template <typename T>
void switchRows ( const int row_src, const int row_des, const int size, T* matrix, T* result )
{
	for (int col = 0; col < size; ++col)
	{
		std::swap(matrix[ index(row_src, col, size) ], matrix[ index(row_des, col, size) ]);
		std::swap(result[ index(row_src, col, size) ], result[ index(row_des, col, size) ]);
	}
}		/* -----  end of template function switchRows  ----- */

template <typename T>
void switchCols ( const int col_src, const int col_des, const int size, T* matrix, T* result )
{
	for (int row = 0; row < size; ++row)
	{
		std::swap(matrix[ index(row, col_src, size) ], matrix[ index(row, col_des, size) ]);
		std::swap(result[ index(row, col_src, size) ], result[ index(row, col_des, size) ]);
	}
}		/* -----  end of template function switchRows  ----- */

template <class T>
void normalizePivotRow ( const int row, const int size, T* matrix, T* result )
{
	T pivot_value = matrix[ index(row, row, size) ];
	for (int col = 0; col < size; ++col)
	{
		matrix[ index(row, col, size) ] /= pivot_value;
		result[ index(row, col, size) ] /= pivot_value;
	}
}		/* -----  end of template function normalizePivotRow  ----- */

template <class T>
void normalizePivotCol ( const int col, const int size, T* matrix, T* result )
{
	T pivot_value = matrix[ index(col, col, size) ];
	for (int row = 0; row < size; ++row)
	{
		matrix[ index(row, col, size) ] /= pivot_value;
		result[ index(row, col, size) ] /= pivot_value;
	}
}		/* -----  end of template function normalizePivotRow  ----- */


template <typename T>
void eliminateByRow ( const int pivot_index, const int size, T* matrix, T* result )
{
	for (int row = 0; row < size; ++row)
	{
		T pivot_col_item = matrix[ index(row, pivot_index, size) ];
		for (int col = 0; col < size; ++col)
		{
			T matrix_pivot_value = matrix[ index(pivot_index, col, size) ];
			T result_pivot_value = result[ index(pivot_index, col, size) ];
			if (row != pivot_index)
			{
				matrix[ index(row, col, size) ] -= pivot_col_item * matrix_pivot_value;
				result[ index(row, col, size) ] -= pivot_col_item * result_pivot_value;
			}
		}
	}
}		/* -----  end of template function eliminateByRow  ----- */

template <typename T>
void eliminateByCol ( const int pivot_index, const int size, T* matrix, T* result )
{
	for (int row = 0; row < size; ++row)
	{
		T matrix_pivot_value = matrix[ index(row, pivot_index, size) ];
		T result_pivot_value = result[ index(row, pivot_index, size) ];
		for (int col = 0; col < size; ++col)
		{
			T pivot_row_item = matrix[ index(pivot_index, col, size) ];
			if (col != pivot_index)
			{
				matrix[ index(row, col, size) ] -= pivot_row_item * matrix_pivot_value;
				result[ index(row, col, size) ] -= pivot_row_item * result_pivot_value;
			}
		}
	}
}		/* -----  end of template function eliminateByRow  ----- */


template <class T>
void genIdentityMatrix ( const int size, T* result )
{
	for (i = 0; i < size; ++i)
	{
		for (j = 0; j < size; ++j)
		{
			if (i == j) {
				result[ index(i, j, size) ] = (T) 1;
			} else {
				result[ index(i, j, size) ] = (T) 0;
			}
		}
	}
}		/* -----  end of template function genIdentityMatrix  ----- */


template <typename T>
void getPivotIndex ( const int index, const int size, T* vector )
{
	int pivot_index = -1;
	int i = index;
	while( pivot_index == -1 && i < size ) {
		if (vector[i] > (T) 0)
		{
			pivot_index = i;
		}
		i++;
	}
	return pivot_index;
}		/* -----  end of template function getPivotIndex  ----- */

template <typename T>
void matrixInvert ( const int size, T* matrix, T* result )
{
	T current_row[size];
	int current_row_size = size * sizeof(T);

	genIdentityMatrix(result, size);
	for (int row = 0; row < size; ++row)
	{
		int idx = row;
		memcpy(&current_row, matrix + row * size, current_row_size);
		int pivot_index = getPivotIndex(idx, size, &current_row);
		if (pivot_index != row)
		{
			switchCols(idx, pivot_index, size, matrix, result);
		}
		normalizePivotRow(row, size, matrix, result);
		eliminateByRow(row, size, matrix, result);
	}
}		/* -----  end of template function matrixInvert  ----- */


template <class T>
void matrixCopyRow ( const int src_row_idx, const int des_row_idx, const int col_num, const T* src, T* des )
{
	for (i = 0; i < col_num; ++i)
	{
		des[ index(des_row_idx, i, col_num) ] = src[ index(src_row_idx, i, col_num) ]; 
	}
}		/* -----  end of template function matrixCopyRow  ----- */
