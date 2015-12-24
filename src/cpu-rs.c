/*
 * =====================================================================================
 *
 *       Filename:  CPU-RS.c
 *
 *    Description:  
 *
 *        Version:  2.0
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
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#define index(i, j, size) (((i) * (size)) + (j))

const int w = 8;

const int NW = 1 << 8;

#define BUFFER_SIZE 256

const unsigned int prim_poly_4 = 023;
const unsigned int prim_poly_8 = 0435;
const unsigned int prim_poly_16 = 0210013;

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
    switch(w) 
    {
        case 4: prim_poly = prim_poly_4; break;
        case 8: prim_poly = prim_poly_8; break;
        case 16: prim_poly = prim_poly_16; break;
        default: return -1;
    }
    x_to_w = 1 << w;
    b = 1;
    r = 0;
    for (log = 0; log < x_to_w-1; log++) 
    {
        if(b>x_to_w) break;
        gflog[b] = (uint8_t) log;
        gfilog[log] = (uint8_t) b;
        b = b << 1;
        if (b & x_to_w) 
        {
            b = b ^ prim_poly;
        }
    }
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
    /* Can't divide by 0 */
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

    for(row=0; row<size; row++)
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
    for(col=0; col<size; col++)
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
    for(row=0; row<size; row++)
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
    for(row=0; row<size; row++)
    {
        pivotColItem = matrix[ index(row, pivotIndex, size) ];
        for(col=0; col<size; col++)
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
    for(row=0; row<size; row++)
    {
        matrixPivotValue = matrix[ index(row, pivotIndex, size) ];
        resultPivotValue = result[ index(row, pivotIndex, size) ];
        for(col=0; col<size; col++)
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
// Gaussian elimination
void invert_matrix(uint8_t *matrix, uint8_t *result, int size)
{
    int row;
    int pivotIndex;
    uint8_t currentRow[size];
    int currentRowSize = size*sizeof(uint8_t);

    get_identity_matrix(result, size);

#ifdef DEBUG
    printf("original matrix:\n");
    show_squre_matrix(matrix, size);
    printf("result:\n");
    show_squre_matrix(result,size);
#endif

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

#ifdef DEBUG
        printf("original matrix:\n");
        show_squre_matrix(matrix, size);
        printf("result:\n");
        show_squre_matrix(result,size);
#endif
        eliminate_by_row(matrix, result, row, size);

#ifdef DEBUG
        printf("original matrix:\n");
        show_squre_matrix(matrix, size);
        printf("result:\n");
        show_squre_matrix(result,size);
#endif
    }
}


void encode_chunk(uint8_t *dataChunk, uint8_t *parityCoeff, uint8_t *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
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

void gen_encoding_matrix(uint8_t *encodingMatrix, int row, int col)
{
    int i;
    int j;
    for(i = 0; i < row; i++)
    {
        for(j = 0; j < col; j++)
        {
            encodingMatrix[i*col + j] = gf_pow(j+1, i);
        }
    }
}

void gen_total_encoding_matrix(uint8_t *totalEncodingMatrix, int nativeBlockNum, int parityBlockNum)
{
    get_identity_matrix(totalEncodingMatrix, nativeBlockNum);
    gen_encoding_matrix(totalEncodingMatrix + nativeBlockNum*nativeBlockNum, parityBlockNum, nativeBlockNum);
}

void write_metadata(char *fileName, int totalSize, int parityBlockNum, int nativeBlockNum)
{
    FILE *fp;
    if( ( fp = fopen(fileName, "wb") ) == NULL )
    {
        printf("Can not open META file!\n");
        exit(0);
    }
    fprintf(fp, "%d\n", totalSize);
    fprintf(fp, "%d %d\n", parityBlockNum, nativeBlockNum);
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
    chunkSize = (totalSize / nativeBlockNum) + (totalSize%nativeBlockNum != 0); 

    uint8_t *dataBuf;
    uint8_t *codeBuf;
    int dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
    int codeSize = parityBlockNum*chunkSize*sizeof(uint8_t);
    dataBuf = (uint8_t*) malloc( nativeBlockNum*chunkSize*sizeof(uint8_t) );
    memset(dataBuf, 0, dataSize);
    codeBuf = (uint8_t*) malloc( parityBlockNum*chunkSize*sizeof(uint8_t) );
    memset(codeBuf, 0, codeSize);

    int i;
    for(i = 0; i < nativeBlockNum; i++)
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
    uint8_t *encodingMatrix;
    encodingMatrix = (uint8_t*) malloc( parityBlockNum*nativeBlockNum*sizeof(uint8_t) );
    gen_encoding_matrix(encodingMatrix, parityBlockNum, nativeBlockNum);
    encode_chunk(dataBuf, encodingMatrix, codeBuf, nativeBlockNum, parityBlockNum, chunkSize);
    clock_gettime(CLOCK_REALTIME,&end);
    totalTime = (double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_nsec-start.tv_nsec)/(double)1000000L;
    printf("Total CPU encoding time: %fms\n", totalTime);
    char metadata_file_name[strlen(file) + 15];
    sprintf(metadata_file_name, "%s.METADATA", file);
    write_metadata(metadata_file_name, totalSize, parityBlockNum, nativeBlockNum);

    char output_file_name[strlen(file) + 5];
    for(i = 0; i < nativeBlockNum; i++)
    {
        //		sprintf(output_file_name, "_%d_", i);
        //		strcat(output_file_name, file);
        sprintf(output_file_name, "_%d_%s", i, file);
        if(( fp_out = fopen(output_file_name, "wb") ) == NULL)
        {
            printf("Can not open source file!\n");
            exit(0);
        }
        if( fwrite(dataBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_out) != sizeof(uint8_t)*chunkSize )
        {
            printf("fwrite error!\n");
            exit(0);
        }
        fclose(fp_out);
    }
    for(i = 0; i < parityBlockNum; i++)
    {
        //		sprintf(output_file_name, "_%d_", i+nativeBlockNum);
        //		strcat(output_file_name, file);
        sprintf(output_file_name, "_%d_%s", i+nativeBlockNum, file);
        if((fp_out = fopen(output_file_name, "wb")) == NULL)
        {
            printf("Can not open source file!\n");
            exit(0);
        }
        if(fwrite(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_out) != sizeof(uint8_t)*chunkSize)
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

// void decode_file(char *inFile, char *confFile, char *outFile, int nativeBlockNum, int parityBlockNum)
void decode_file(char *inFile, char *confFile, char *outFile)
{
    int chunkSize = 1;
    int totalSize;

    int parityBlockNum;
    int nativeBlockNum;

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
    char metadata_file_name[strlen(inFile) + 15];
    sprintf(metadata_file_name, "%s.METADATA", inFile);
    if((fp_in = fopen(metadata_file_name, "rb")) == NULL)
    {
        printf("Can not open metadata file!\n");
        exit(0);
    }
    fscanf(fp_in, "%d", &totalSize);
    fscanf(fp_in, "%d %d", &parityBlockNum, &nativeBlockNum);
    fclose(fp_in);
    //	chunkSize = (int) (ceil( (float) (totalSize / nativeBlockNum) )); 
    chunkSize = (totalSize / nativeBlockNum) + ( totalSize%nativeBlockNum != 0 ); 

#ifdef DEBUG
    printf("chunk size: %d\n", chunkSize);
#endif

    totalMatrixSize = nativeBlockNum * (nativeBlockNum + parityBlockNum);
    totalEncodingMatrix = (uint8_t*) malloc(totalMatrixSize);
    matrixSize = nativeBlockNum * nativeBlockNum;
    encodingMatrix = (uint8_t*) malloc(matrixSize);
    gen_total_encoding_matrix(totalEncodingMatrix, nativeBlockNum, parityBlockNum);

    dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
    codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
    dataBuf = (uint8_t*) malloc( dataSize );
    memset(dataBuf, 0, dataSize);
    codeBuf = (uint8_t*) malloc( codeSize);
    memset(codeBuf, 0, codeSize);

    FILE *fp_conf;
    char input_file_name[strlen(inFile) + 20];
    int index;
    fp_conf = fopen(confFile, "r");
    int i;
    for(i = 0; i < nativeBlockNum; i++)
    {
        fscanf(fp_conf, "%s", input_file_name);
        index = atoi(input_file_name + 1);

        copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

        fp_in = fopen(input_file_name, "rb");
        fseek(fp_in, 0L, SEEK_SET);
        // this part can be process in parallel with computing inversed matrix
        fread(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
        fclose(fp_in);
    }
    fclose(fp_conf);

    struct timespec start, end;
    double totalTime;
    clock_gettime(CLOCK_REALTIME,&start);
    uint8_t *decodingMatrix;
    decodingMatrix = (uint8_t*) malloc( matrixSize );

    invert_matrix(encodingMatrix, decodingMatrix, nativeBlockNum);

    //#ifndef DEBUG
    //	show_matrix(totalEncodingMatrix, nativeBlockNum+parityBlockNum, nativeBlockNum);
    //#endif

    decode_chunk(dataBuf, decodingMatrix, codeBuf, nativeBlockNum, parityBlockNum, chunkSize);
    clock_gettime(CLOCK_REALTIME,&end);
    totalTime = (double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_nsec-start.tv_nsec)/(double)1000000L;
    printf("Total CPU decoding time: %fms\n", totalTime);

    if(outFile == NULL)
    {
        fp_out = fopen(inFile, "wb");
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

void show_help_info()
{
    printf("Usage:\n");
    printf("[-h]: show usage information\n");
    printf("Encode: [-k|-K nativeBlockNum] [-n|-N totalBlockNum] [-e|-E fileName]\n");
    printf("Decode: [-d|-D] [-k|-K nativeBlockNum] [-n|-N totalBlockNum] \n\t [-i|-I originalFileName] [-c|-C config] [-o|-O output]\n");
    printf("For encoding, the -k, -n, and -e options are all necessary.\n");
    printf("For decoding, the -d, -i, and -c options are all necessary.\n");
    printf("If the -o option is not set, the original file name will be chosen as the output file name by default.\n");
    exit(0);
}

int main(int argc, char *argv[])
{
    int nativeBlockNum = 0;
    int parityBlockNum = 0;
    int totalBlockNum = 0;
    char *inFile = NULL;
    char *confFile = NULL;
    char *outFile = NULL;

    enum func
    {
        encode,
        decode
    };
    enum func op;
    int func_flag = 0;

    int option;
    while((option = getopt(argc, argv, "Kk:Nn:Ee:Ii:Cc:Oo:Ddh")) != -1) {
        switch ( option ) {
            case 'K':	
            case 'k':	
                nativeBlockNum = (int) atoi(optarg);
                break;

            case 'N':	
            case 'n':	
                totalBlockNum = (int) atoi(optarg);
                break;

            case 'E':	
            case 'e':	
                inFile = optarg;
                op = encode;
                func_flag = 1;
                break;

            case 'D':	
            case 'd':	
                op = decode;
                func_flag = 1;
                break;

            case 'I':	
            case 'i':	
                if (func_flag == 1 && op == decode)
                {
                    inFile = optarg;
                }
                else
                {
                    show_help_info();
                }
                break;

            case 'C':	
            case 'c':	
                if (func_flag == 1 && op == decode)
                {
                    confFile = optarg;
                }
                else
                {
                    show_help_info();
                }
                break;

            case 'O':	
            case 'o':	
                if (func_flag == 1 && op == decode)
                {
                    outFile = optarg;
                }
                else
                {
                    show_help_info();
                }
                break;

            case 'h':	
                show_help_info();
                break;

            default:	
                show_help_info();
                break;
        }	/* -----  end switch  ----- */
    }

    // setup table for GF(2^8)
    setup_tables(8);

    switch ( op ) {
        case encode:	
            assert(nativeBlockNum != 0);
            assert(totalBlockNum != 0);
            parityBlockNum = totalBlockNum - nativeBlockNum;
            encode_file(inFile, nativeBlockNum, parityBlockNum);
            break;

        case decode:	
            assert(inFile != NULL);
            assert(confFile != NULL);
            assert(outFile != NULL);
            decode_file(inFile, confFile, outFile);
            break;

        default:	
            break;
    }		/* -----  end switch  ----- */

    return 0;

}
