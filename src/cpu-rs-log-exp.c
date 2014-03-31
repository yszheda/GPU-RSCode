/*
 * =====================================================================================
 *
 *       Filename:  CPU-RS.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/05/2012 10:42:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuai YUAN(yszheda AT gmail.com), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include "galoisfield.h"

#define IDC2D(i,j,ld) (((i)*(ld))+(j))

//extern uint8_t encodingMatrix[];
//extern unsigned char encodingMatrix[];

//unsigned char* encodingMatrix;

const int w = 8;

const int NW = 1 << 8;
//#define NW (1 << w) /* In other words, NW equals 2 to the w-th power */

#define DEBUG 
#define BUFFER_SIZE 256

//unsigned int prim_poly_4 = 023;
//unsigned int prim_poly_8 = 0435;
//unsigned int prim_poly_16 = 0210013;
unsigned int prim_poly_4 = 023;
unsigned int prim_poly_8 = 0435;
//uint8_t prim_poly_8 = 285;
unsigned int prim_poly_16 = 0210013;
//uint8_t *gflog;
//uint8_t *gfilog;
uint8_t gflog[256];
uint8_t gfilog[256];

int setup_tables(int w)
{
	unsigned int b;
	unsigned int r;
   	unsigned int log;
	unsigned int x_to_w;
	unsigned int prim_poly;
	unsigned int x;
	unsigned int y;
//	uint8_t b, log, x_to_w, prim_poly;
	switch(w) 
	{
		case 4: prim_poly = prim_poly_4; break;
		case 8: prim_poly = prim_poly_8; break;
		case 16: prim_poly = prim_poly_16; break;
		default: return -1;
	}
	x_to_w = 1 << w;
//	gflog = (uint8_t *) malloc (sizeof(uint8_t) * x_to_w);
//	gfilog = (uint8_t *) malloc (sizeof(uint8_t) * x_to_w);
	b = 1;
	r = 0;
	for (log = 0; log < x_to_w-1; log++) 
	{
		/*
		r = 0;
		x = 1;
		y = log;
		while(y)
		{
			printf("y=%d\n",y);
			if(y & 1)
			{
				r = r ^ b;
			}
			y = y >> 1;
			x = x << 1;
			if (x & x_to_w) x = x ^ prim_poly;
		}
			printf("r=%d\n",r);
			printf("log=%d\n",log);
		*/
		if(b>x_to_w) break;
		gflog[b] = (uint8_t) log;
		gfilog[log] = (uint8_t) b;
		b = b << 1;
		if (b & x_to_w) 
		{
			b = b ^ prim_poly;
//#ifdef DEBUG
//printf("prim_poly=%d\n", prim_poly);
//printf("test b=%d\n", b);
//#endif
		}
	}
//#ifdef DEBUG
//	printf("b=%d\n",b);
//#endif
	return 0;
}

void showgflog()
{
	int i;
	printf("gflog\n");
	for(i =0; i< NW; i++)
	{
		printf("%d ", gflog[i]);
	}
}

void showgfilog()
{
	int i;
	printf("gfilog\n");
	for(i =0; i< NW; i++)
	{
		printf("%d ", gfilog[i]);
	}
}

uint8_t gf_add(uint8_t a, uint8_t b)
{
	return a^b;
}

uint8_t gf_sub(uint8_t a, uint8_t b)
{
	return gf_add(a, b);
}

uint8_t gf_mul(uint8_t a, uint8_t b)
{
	int sum_log;
	if (a == 0 || b == 0)
	{
		return 0;
	}
//	sum_log = (gflog[a] + gflog[b]) % (NW-1);
	sum_log = gflog[a] + gflog[b];
	if (sum_log >= NW-1)
	{	
		sum_log -= NW-1;
	}
	return gfilog[sum_log];
}

uint8_t gf_mul_bit(uint8_t a, uint8_t b)
{
	int sum_log;
	while(b)
	{
		if(b & 1)
		{
			sum_log ^= a;
		}
		a = (a << 1) ^ (a & 0x80? 0x1d: 0);
		b >>= 1;
	}
	return sum_log;
}

uint8_t gf_div(uint8_t a, uint8_t b)
{
	int diff_log;
	if (a == 0)
	{	
		return 0;
	}
	/* Can’t divide by 0 */
	if (b == 0)
	{
		return -1;
	}
//	diff_log = (gflog[a] - gflog[b]) % (NW-1);
	diff_log = gflog[a] - gflog[b];
	if (diff_log < 0)
	{	
		diff_log += NW-1;
	}
	return gfilog[diff_log];
}

uint8_t gf_pow(uint8_t a, uint8_t power)
{
	int pow_log = (gflog[a] * power) % (NW-1);
	return gfilog[pow_log];
}

void show_squre_matrix(uint8_t *matrix, int size)
{
	int i;
	int j;
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			printf("%d ", matrix[i*size+j]);
		}
		printf("\n");
	}
}

void gen_encoding_matrix(uint8_t *encodingMatrix, int row, int col)
{
	int i;
	int j;
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			encodingMatrix[i*col + j] = gf_pow((j+1) % 256, i);
		}
	}
}

// C=AB
// A: nxp
// B: pxm
// C: nxm
void matrix_mul(uint8_t *A, uint8_t *B, uint8_t *C, int n, int p, int m)
{
	int i;
	int j;
	int k;
	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
		{
			for(k=0; k<p; k++)
			{
//				C[i*m+j] = gf_add(C[i*m+j], gf_mul(A[i*p+k],B[k*m+j]));
				C[i*m+j] ^= gf_mul(A[i*p+k],B[k*m+j]);
			}
		}
	}
}

// switch rows if the current row is not the pivot row
void switch_rows(uint8_t *matrix, uint8_t *result, int rowSrc, int rowDes, int size)
{
	int col;
    uint8_t oldMatrixItem;
    uint8_t oldResultItem;

	for(col=0; col<size; col++)
    {
        oldMatrixItem = matrix[ rowSrc*size+col ];
        matrix[ rowSrc*size+col ] = matrix[ rowDes*size+col ];
        matrix[ IDC2D(rowDes, col, size) ] = oldMatrixItem; 

        oldResultItem = result[ IDC2D(rowSrc, col, size) ];
        result[ IDC2D(rowSrc, col, size) ] = result[ IDC2D(rowDes, col, size) ];
        result[ IDC2D(rowDes, col, size) ] = oldResultItem; 
    }
} 

void switch_columns(uint8_t *matrix, uint8_t *result, int colSrc, int colDes, int size)
{
    int row;
    uint8_t oldMatrixItem;
    uint8_t oldResultItem;

	for(row=0; row<size; row++)
    {
        oldMatrixItem = matrix[ IDC2D(row, colSrc, size) ];
        matrix[ IDC2D(row, colSrc, size) ] = matrix[ IDC2D(row, colDes, size) ];
        matrix[ IDC2D(row, colDes, size) ] = oldMatrixItem; 

        oldResultItem = result[ IDC2D(row, colSrc, size) ];
        result[ IDC2D(row, colSrc, size) ] = result[ IDC2D(row, colDes, size) ];
        result[ IDC2D(row, colSrc, size) ] = oldResultItem; 
    }
} 

// normalize the row by the pivot value
void normalize_pivot_row(uint8_t *matrix, uint8_t *result, int row, int size)
{
	int col;
    uint8_t pivotValue;

    pivotValue = matrix[ IDC2D(row, row, size) ];
	for(col=0; col<size; col++)
    {
        matrix[ IDC2D(row, col, size)] = gf_div(matrix[ IDC2D(row, col, size) ], pivotValue);
        result[ IDC2D(row, col, size)] = gf_div(result[ IDC2D(row, col, size) ], pivotValue);
    }
}

// normalize the column by the pivot value
void normalize_pivot_col(uint8_t *matrix, uint8_t *result, int col, int size)
{
	int row;
    uint8_t pivotValue;

    pivotValue = matrix[ IDC2D(col, col, size) ];
	for(row=0; row<size; row++)
    {
        matrix[ IDC2D(row, col, size)] = gf_div(matrix[ IDC2D(row, col, size) ], pivotValue);
        result[ IDC2D(row, col, size)] = gf_div(result[ IDC2D(row, col, size) ], pivotValue);
    }
}

//eliminate by row to make the pivot column become reduced echelon form
void eliminate_by_row(uint8_t *matrix, uint8_t *result, int pivotIndex, int size)
{
	int row;
	int col;
	uint8_t matrixPivotValue;
	uint8_t resultPivotValue;
	uint8_t pivotColItem;
	for(row=0; row<size; row++)
	{
		pivotColItem = matrix[ IDC2D(row, pivotIndex, size) ];
		for(col=0; col<size; col++)
		{
            matrixPivotValue = matrix[ IDC2D(pivotIndex, col, size) ];
            resultPivotValue = result[ IDC2D(pivotIndex, col, size) ];
			if(row != pivotIndex)
			{
				matrix[ IDC2D(row, col, size) ] ^= gf_mul(pivotColItem, matrixPivotValue);
				result[ IDC2D(row, col, size) ] ^= gf_mul(pivotColItem, resultPivotValue);
			}
		}
	}
}

//eliminate by column to make the pivot row become reduced echelon form
void eliminate_by_col(uint8_t *matrix, uint8_t *result, int pivotIndex, int size)
{
	int row;
	int col;
	uint8_t matrixPivotValue;
	uint8_t resultPivotValue;
	uint8_t pivotRowItem;
	for(row=0; row<size; row++)
	{
        matrixPivotValue = matrix[ IDC2D(row, pivotIndex, size) ];
        resultPivotValue = result[ IDC2D(row, pivotIndex, size) ];
		for(col=0; col<size; col++)
		{
			pivotRowItem = matrix[ IDC2D(pivotIndex, col, size) ];
			if(col != pivotIndex)
			{
				matrix[ IDC2D(row, col, size) ] ^= gf_mul(pivotRowItem, matrixPivotValue);
				result[ IDC2D(row, col, size) ] ^= gf_mul(pivotRowItem, resultPivotValue);
			}
		}
	}
}


//generate an identity matrix
void get_identity_matrix(uint8_t *result, int size)
{
	int i;
	int j;
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			if(i == j)
			{
				result[i*size+j] = 1;
			}
			else
			{
				result[i*size+j] = 0;
			}
		}
	}
}


//find the pivot index in the given row/column
int get_pivot_index(uint8_t *vector, int index, int size)
{
    int pivotIndex = -1;
    int i = index;
    while( pivotIndex == -1 && i < size )
    {
        pivotIndex = (vector[i] > 0)? i: -1;        
        i++;
    }
    return pivotIndex;
}

// compute the inverse of a given matrix
// Guassian elimination
void invert_matrix(uint8_t *matrix, uint8_t *result, int size)
{
	int row;
	int pivotIndex;
    uint8_t currentRow[size];
    int currentRowSize = size*sizeof(uint8_t);

    get_identity_matrix(result, size);
	
// #ifdef DEBUG
// printf("original matrix:\n");
// 	show_squre_matrix(matrix, size);
// printf("result: \n");
// show_squre_matrix(result,size);
// #endif
	for(row=0; row<size; row++)
    {
		// check whether the leading coefficient of the current row is in the 'index'th column
		int index = row;
        memcpy(&currentRow, matrix+row*size, currentRowSize);
        pivotIndex = get_pivot_index(currentRow, index, size);
        if( pivotIndex != row )
		{
            switch_columns(matrix, result, index, pivotIndex, size);
		}

    	// Normalize the pivot row
        normalize_pivot_row(matrix, result, index, size);

// #ifdef DEBUG
// printf("original matrix:\n");
// 	show_squre_matrix(matrix, size);
// printf("result: \n");
// show_squre_matrix(result,size);
// #endif
        eliminate_by_row(matrix, result, row, size);

// #ifdef DEBUG
// printf("original matrix:\n");
// 	show_squre_matrix(matrix, size);
// printf("result: \n");
// show_squre_matrix(result,size);
// #endif
    }
#ifdef DEBUG
// printf("original matrix:\n");
//	show_squre_matrix(matrix, size);
printf("result: \n");
show_squre_matrix(result,size);
#endif
}


void encode_chunk(uint8_t *dataChunk, uint8_t *parityCoeff, uint8_t *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
//	codeChunk = (unsigned char*)malloc(parityBlockNum*chunkSize);
//	codeChunk = (uint8_t*) malloc( parityBlockNum*chunkSize*sizeof(uint8_t) );
	matrix_mul(parityCoeff, dataChunk, codeChunk, parityBlockNum, nativeBlockNum, chunkSize);
}

void decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	matrix_mul(parityCoeff, codeChunk, dataChunk, nativeBlockNum, nativeBlockNum, chunkSize);
}

void show_code_chunk(unsigned char *codeChunk, int parityBlockNum, int chunkSize)
{
	int i;
	int j;
	for(i=0; i<parityBlockNum; i++)
	{
		for(j=0; j<chunkSize; j++)
		{
			printf("%d ", codeChunk[i*chunkSize+j]);
		}
		printf("\n");
	}
}

void show_matrix(uint8_t *matrix, int row, int col)
{
	int i;
	int j;
	for(i=0; i<row; i++)
	{
		for(j=0; j<col; j++)
		{
			printf("%d ", matrix[i*col+j]);
		}
		printf("\n");
	}
}

void copy_matrix(uint8_t *src, uint8_t *des, int srcRowIndex, int desRowIndex, int rowSize)
{
	int i;
	for(i=0; i<rowSize; i++)
	{
		des[desRowIndex*rowSize+i] = src[srcRowIndex*rowSize+i];
	}
}

void write_metadata(int totalSize, int parityBlockNum, int nativeBlockNum, uint8_t *encodingMatrix)
{
	FILE *fp;
	if( ( fp = fopen(".METADATA", "wb") ) == NULL )
	{
		printf("Can not open META file!\n");
		exit(0);
	}
	fprintf(fp, "%d\n", totalSize);
	fprintf(fp, "%d %d\n", parityBlockNum, nativeBlockNum);
	int i;
	int j;
	for(i=0; i<nativeBlockNum; i++)
	{
		for(j=0; j<nativeBlockNum; j++)
		{
			if(i == j)
			{
				fprintf(fp, "1 ");
			}
			else
			{
				fprintf(fp, "0 ");
			}
		}
		fprintf(fp, "\n");
	}
	for(i=0; i<parityBlockNum; i++)
	{
		for(j=0; j<nativeBlockNum; j++)
		{
			fprintf(fp, "%d ", encodingMatrix[i*nativeBlockNum+j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}


void encode_file(char *file, int nativeBlockNum, int parityBlockNum)
{
	int chunkSize = 1;
	int totalSize;

	FILE *fp_in;
	FILE *fp_out;
	if( ( fp_in = fopen(file,"rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}

	fseek(fp_in, 0L, SEEK_END);
	//ftell() get the total size of the file
	totalSize = ftell(fp_in);
	chunkSize = (totalSize / nativeBlockNum) + ( totalSize%nativeBlockNum != 0 ); 
//	chunkSize = (ftell(fp_in) / nativeBlockNum) + ( ftell(fp_in)%nativeBlockNum != 0 ); 
//	chunkSize = (int) (ceil( (long double) (ftell(fp_in) / nativeBlockNum)) ); 

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host
	int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	int codeSize = parityBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( nativeBlockNum*chunkSize*sizeof(uint8_t) );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( parityBlockNum*chunkSize*sizeof(uint8_t) );
	memset(codeBuf, 0, codeSize);
	
	int i;
	for(i=0; i<nativeBlockNum; i++)
	{
		if( fseek( fp_in, i*chunkSize, SEEK_SET ) == -1 )
		{
			printf("fseek error!\n");
			exit(0);
		}

		if( fread( dataBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in ) == EOF )
		{
			printf("fread error!\n");
			exit(0);
		}
	}
	fclose(fp_in);
	
	struct timespec start, end;
	double totalTime;
	clock_gettime(CLOCK_REALTIME,&start);
//	// setup table for GF(2^8)
//	setup_tables(8);
	uint8_t *encodingMatrix;
	encodingMatrix = (uint8_t*) malloc( parityBlockNum*nativeBlockNum*sizeof(uint8_t) );
	gen_encoding_matrix(encodingMatrix, parityBlockNum, nativeBlockNum);
	write_metadata(totalSize, parityBlockNum, nativeBlockNum, encodingMatrix);
	encode_chunk(dataBuf, encodingMatrix, codeBuf, nativeBlockNum, parityBlockNum, chunkSize);
//	matrix_mul(encodingMatrix, dataBuf, codeBuf, parityBlockNum, nativeBlockNum, chunkSize);
/*
	//int i;
	int j;
	int k;
	int n=parityBlockNum;
	int p=nativeBlockNum;
	int m=chunkSize;
	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
		{
			codeBuf[i*m+j] = 0;
			for(k=0; k<p; k++)
			{
//				C[i*m+j] = gf_add(C[i*m+j], gf_mul(A[i*p+k],B[k*m+j]));
				codeBuf[i*m+j] ^= gf_mul(encodingMatrix[i*p+k],dataBuf[k*m+j]);
			}
		}
	}
*/
	clock_gettime(CLOCK_REALTIME,&end);
	totalTime = (double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_nsec-start.tv_nsec)/(double)1000000L;
	printf("Total CPU encoding time: %fms\n", totalTime);

	char output_file_name[100];
	for(i=0; i<nativeBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_", i);
		strcat(output_file_name, file);
		if( ( fp_out = fopen(output_file_name, "wb") ) == NULL )
		{
			printf("Can not open source file!\n");
			exit(0);
		}
		if( fwrite(dataBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_out ) != sizeof(uint8_t)*chunkSize )
		{
			printf("fwrite error!\n");
			exit(0);
		}
		fclose(fp_out);
	}
	for(i=0; i<parityBlockNum; i++)
	{
		sprintf(output_file_name, "_%d_", i+nativeBlockNum);
		strcat(output_file_name, file);
		if( ( fp_out = fopen(output_file_name, "wb") ) == NULL )
		{
			printf("Can not open source file!\n");
			exit(0);
		}
		if( fwrite(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_out ) != sizeof(uint8_t)*chunkSize )
		{
			printf("fwrite error!\n");
			exit(0);
		}
		fclose(fp_out);
	}

	free(dataBuf);
	free(codeBuf);

	free(encodingMatrix);

}

void decode_file(char *confFile, char *outFile, int nativeBlockNum, int parityBlockNum)
{
	int chunkSize = 1;
	int totalSize;

	uint8_t *dataBuf;
	uint8_t *codeBuf;

	int dataSize;
	int codeSize;

	FILE *fp_in;
	FILE *fp_out;

	int totalMatrixSize;
	int matrixSize;
	uint8_t *totalEncodingMatrix;
	uint8_t *encodingMatrix;
	if( ( fp_in = fopen(".METADATA","rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}
	fscanf(fp_in, "%d", &totalSize);
	fscanf(fp_in, "%d %d", &parityBlockNum, &nativeBlockNum);
//	chunkSize = (int) (ceil( (float) (totalSize / nativeBlockNum) )); 
	chunkSize = (totalSize / nativeBlockNum) + ( totalSize%nativeBlockNum != 0 ); 
#ifdef DEBUG
printf("chunk size: %d\n", chunkSize);
#endif
	totalMatrixSize = nativeBlockNum * ( nativeBlockNum + parityBlockNum );
	totalEncodingMatrix = (uint8_t*) malloc( totalMatrixSize );
	matrixSize = nativeBlockNum * nativeBlockNum;
	encodingMatrix = (uint8_t*) malloc( matrixSize );
	int i;
	for(i =0; i<nativeBlockNum*(nativeBlockNum+parityBlockNum); i++)
	{
		fscanf(fp_in, "%d", totalEncodingMatrix+i);
	}

	dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( dataSize );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( codeSize);
	memset(codeBuf, 0, codeSize);

	if(confFile != NULL)
	{
		FILE *fp_conf;
		char input_file_name[100];
		int index;
		fp_conf = fopen(confFile, "r");

		for(i=0; i<nativeBlockNum; i++)
		{
			fscanf(fp_conf, "%s", input_file_name);
			index = atoi(input_file_name+1);

			copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

			fp_in = fopen(input_file_name, "rb");
			fseek(fp_in, 0L, SEEK_SET);
			// this part can be process in parallel with computing inversed matrix
			fread(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
			fclose(fp_in);
		}
		fclose(fp_conf);
	}
	else
	{
		for(i=0; i<nativeBlockNum; i++)
		{
			char input_file_name[100];
			int index;
			printf("Please enter the file name of fragment:\n");
			scanf("%s", input_file_name);
			index = atoi(input_file_name+1);
			printf("#%dth fragment\n", index);

			copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

			fp_in = fopen(input_file_name, "rb");
			fseek(fp_in, 0L, SEEK_SET);
			// TODO: this part can be process in parallel with computing inversed matrix
			fread(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
			fclose(fp_in);

		}
	}
	
	struct timespec start, end;
	double totalTime;
	clock_gettime(CLOCK_REALTIME,&start);
	uint8_t *decodingMatrix;
	decodingMatrix = (uint8_t*) malloc( matrixSize );

    invert_matrix(encodingMatrix, decodingMatrix, nativeBlockNum);
//#ifndef DEBUG
//	show_matrix(totalEncodingMatrix, nativeBlockNum+parityBlockNum, nativeBlockNum);
//#endif

//#ifndef DEBUG
	decode_chunk(dataBuf, decodingMatrix, codeBuf, nativeBlockNum, parityBlockNum, chunkSize);
//#endif
//#ifdef DEBUG
//	uint8_t test_DM[16] = {1,0,0,0, 2,1,3,7, 3,1,2,6, 0,0,0,1};	
//	decode_chunk(dataBuf, test_DM, codeBuf, nativeBlockNum, parityBlockNum, chunkSize);
//#endif	
	clock_gettime(CLOCK_REALTIME,&end);
	totalTime = (double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_nsec-start.tv_nsec)/(double)1000000L;
	printf("Total CPU decoding time: %fms\n", totalTime);

	if(outFile == NULL)
	{
		char output_file_name[100];
		printf("Enter the name of the decoded file:\n");
		scanf("%s", output_file_name);
		fp_out = fopen(output_file_name, "wb");
	}
	else
	{
		fp_out = fopen(outFile, "wb");
	}
	fwrite(dataBuf, sizeof(uint8_t), totalSize, fp_out);
	fclose(fp_out);

	free(dataBuf);
	free(codeBuf);

}

int main(int argc, char *argv[])
{
	int nativeBlockNum = 4;
	int parityBlockNum = 2;
	char *inFile = NULL;
	char *confFile = NULL;
	char *outFile = NULL;

	enum func
	{
		encode,
		decode
	};
	enum func op;

	nativeBlockNum = atoi(argv[1]);
	parityBlockNum = atoi(argv[2]);

	if( strcmp(argv[3], "-e") == 0 )
	{
		op = encode;
	}
	else if( strcmp(argv[3], "-d") == 0 )
	{
		op = decode;
	}
	else
	{
		printf("Invalid option!\n");
		exit(-1);
	}
	
	// setup table for GF(2^8)
	setup_tables(8);

	switch(op)
	{
		case encode:
			inFile = argv[4];
			encode_file(inFile, nativeBlockNum, parityBlockNum);
			break;

		case decode:
			if(argc == 5)
			{
				confFile = argv[4];
			}
			else if(argc == 7 && strcmp(argv[5], "-o") == 0)
			{
				confFile = argv[4];
				outFile = argv[6];
			}
			else
			{
				printf("Invalid command!\n");
				exit(-1);
			}
			decode_file(confFile, outFile, nativeBlockNum, parityBlockNum);
			break;
	}

	return 0;

}
