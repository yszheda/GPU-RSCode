/*
 * =====================================================================================
 *
 *       Filename:  cpu-decode.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/04/14 08:49:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuai YUAN(yszheda AT gmail.com), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "cpu-decode.h"

const int width = 8;
const int field_size = 1 << 8;

uint8_t gflog[256];
uint8_t gfexp[256];

int setup_host_tables()
{
	const int width = 8;
	const int field_size = 1 << width;
	const unsigned int prim_poly = 0435;
   	int log;
	int exp = 1;
	// use int as book-keeping index instead of unsigned int
	for (log = 0; log < field_size - 1; log++) 
	{
		if(exp > field_size) break;
		gflog[exp] = (uint8_t) log;
		gfexp[log] = (uint8_t) exp;
		exp = exp << 1;
		if (exp & field_size) 
		{
			exp = exp ^ prim_poly;
		}
	}
	return 0;
}

inline uint8_t gf_add(uint8_t a, uint8_t b)
{
	return a ^ b;
}

inline uint8_t gf_sub(uint8_t a, uint8_t b)
{
	return gf_add(a, b);
}

inline uint8_t gf_mul(uint8_t a, uint8_t b)
{
	int sum_log;
	if (a == 0 || b == 0)
	{
		return 0;
	}
//	sum_log = (gflog[a] + gflog[b]) % (field_size-1);
	sum_log = gflog[a] + gflog[b];
	if (sum_log >= field_size - 1)
	{	
		sum_log -= field_size - 1;
	}
	return gfexp[sum_log];
}

inline uint8_t gf_div(uint8_t a, uint8_t b)
{
	int diff_log;
	if (a == 0)
	{	
		return 0;
	}
	// optimize out exception cases
	/*
	// Can't divide by 0
	if (b == 0)
	{
		return -1;
	}
	*/
//	diff_log = (gflog[a] - gflog[b]) % (field_size-1);
	diff_log = gflog[a] - gflog[b];
	if (diff_log < 0)
	{	
		diff_log += field_size - 1;
	}
	return gfexp[diff_log];
}

// switch rows if the current row is not the pivot row
void switch_rows(uint8_t *matrix, uint8_t *result, int rowSrc, int rowDes, int size)
{
	int col;
    uint8_t oldMatrixItem;
    uint8_t oldResultItem;

	for(col = 0; col < size; col++)
    {
        oldMatrixItem = matrix[ index(rowSrc, col, size) ];
        matrix[ index(rowSrc, col, size) ] = matrix[ index(rowDes, col, size) ];
        matrix[ index(rowDes, col, size) ] = oldMatrixItem; 

        oldResultItem = result[ index(rowSrc, col, size) ];
        result[ index(rowSrc, col, size) ] = result[ index(rowDes, col, size) ];
        result[ index(rowDes, col, size) ] = oldResultItem; 
    }
} 

void switch_columns(uint8_t *matrix, uint8_t *result, int colSrc, int colDes, int size)
{
    int row;
    uint8_t oldMatrixItem;
    uint8_t oldResultItem;

	for(row = 0; row < size; row++)
    {
        oldMatrixItem = matrix[ index(row, colSrc, size) ];
        matrix[ index(row, colSrc, size) ] = matrix[ index(row, colDes, size) ];
        matrix[ index(row, colDes, size) ] = oldMatrixItem; 

        oldResultItem = result[ index(row, colSrc, size) ];
        result[ index(row, colSrc, size) ] = result[ index(row, colDes, size) ];
        result[ index(row, colSrc, size) ] = oldResultItem; 
    }
} 

// normalize the row by the pivot value
void normalize_pivot_row(uint8_t *matrix, uint8_t *result, int row, int size)
{
	int col;
    uint8_t pivotValue;

    pivotValue = matrix[ index(row, row, size) ];
	for(col = 0; col < size; col++)
    {
        matrix[ index(row, col, size)] = gf_div(matrix[ index(row, col, size) ], pivotValue);
        result[ index(row, col, size)] = gf_div(result[ index(row, col, size) ], pivotValue);
    }
}

// normalize the column by the pivot value
void normalize_pivot_col(uint8_t *matrix, uint8_t *result, int col, int size)
{
	int row;
    uint8_t pivotValue;

    pivotValue = matrix[ index(col, col, size) ];
	for(row = 0; row < size; row++)
    {
        matrix[ index(row, col, size)] = gf_div(matrix[ index(row, col, size) ], pivotValue);
        result[ index(row, col, size)] = gf_div(result[ index(row, col, size) ], pivotValue);
    }
}

// eliminate by row to make the pivot column become reduced echelon form
void eliminate_by_row(uint8_t *matrix, uint8_t *result, int pivotIndex, int size)
{
	int row;
	int col;
	uint8_t matrixPivotValue;
	uint8_t resultPivotValue;
	uint8_t pivotColItem;
	for(row = 0; row < size; row++)
	{
		pivotColItem = matrix[ index(row, pivotIndex, size) ];
		for(col = 0; col < size; col++)
		{
            matrixPivotValue = matrix[ index(pivotIndex, col, size) ];
            resultPivotValue = result[ index(pivotIndex, col, size) ];
			if(row != pivotIndex)
			{
				matrix[ index(row, col, size) ] ^= gf_mul(pivotColItem, matrixPivotValue);
				result[ index(row, col, size) ] ^= gf_mul(pivotColItem, resultPivotValue);
			}
		}
	}
}

// eliminate by column to make the pivot row become reduced echelon form
void eliminate_by_col(uint8_t *matrix, uint8_t *result, int pivotIndex, int size)
{
	int row;
	int col;
	uint8_t matrixPivotValue;
	uint8_t resultPivotValue;
	uint8_t pivotRowItem;
	for(row = 0; row < size; row++)
	{
        matrixPivotValue = matrix[ index(row, pivotIndex, size) ];
        resultPivotValue = result[ index(row, pivotIndex, size) ];
		for(col = 0; col < size; col++)
		{
			pivotRowItem = matrix[ index(pivotIndex, col, size) ];
			if(col != pivotIndex)
			{
				matrix[ index(row, col, size) ] ^= gf_mul(pivotRowItem, matrixPivotValue);
				result[ index(row, col, size) ] ^= gf_mul(pivotRowItem, resultPivotValue);
			}
		}
	}
}

// generate an identity matrix
void get_identity_matrix(uint8_t *result, int size)
{
	int i;
	int j;
	for(i = 0; i < size; i++)
	{
		for(j = 0; j < size; j++)
		{
			if(i == j)
			{
				result[ index(i, j, size) ] = 1;
			}
			else
			{
				result[ index(i, j, size) ] = 0;
			}
		}
	}
}

//find the pivot index in the given row/column
int get_pivot_index(uint8_t *vector, int index, int size)
{
    int pivotIndex = -1;
    int i = index;
    while(pivotIndex == -1 && i < size)
    {
        pivotIndex = (vector[i] > 0)? i: -1;        
        i++;
    }
    return pivotIndex;
}

// compute the inverse of a given matrix
// Gaussian elimination
void CPU_invert_matrix(uint8_t *matrix, uint8_t *result, int size)
{
	int row;
	int pivotIndex;
    uint8_t currentRow[size];
    int currentRowSize = size * sizeof(uint8_t);
	
	setup_host_tables();

    get_identity_matrix(result, size);
	
#ifdef DEBUG
printf("original matrix:\n");
show_squre_matrix(matrix, size);
printf("result:\n");
show_squre_matrix(result, size);
#endif

	for(row = 0; row < size; row++)
    {
		// check whether the leading coefficient of the current row is in the 'index'th column
		int index = row;
        memcpy(&currentRow, matrix + row * size, currentRowSize);
        pivotIndex = get_pivot_index(currentRow, index, size);
        if(pivotIndex != row)
		{
            switch_rows(matrix, result, index, pivotIndex, size);
		}

    	// Normalize the pivot row
        normalize_pivot_row(matrix, result, index, size);

#ifdef DEBUG
printf("original matrix:\n");
show_squre_matrix(matrix, size);
printf("result:\n");
show_squre_matrix(result, size);
#endif
        eliminate_by_row(matrix, result, row, size);

#ifdef DEBUG
printf("original matrix:\n");
show_squre_matrix(matrix, size);
printf("result:\n");
show_squre_matrix(result, size);
#endif
    }
}

