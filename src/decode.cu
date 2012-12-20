#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <stdint.h>
//#include "galoisfield.h"

#define SQUARE_BLOCK_SIZE 16    // MAX 
#define LINEAR_BLOCK_SIZE 512   // MAX 

#define DISPLAY_SETTINGS false
#define DISPLAY false

//#define IDC2D(i,j,ld) (((j)*(ld))+(i))
#define IDC2D(i,j,ld) (((i)*(ld))+(j))

#define W 8
#define NW (1 << W) /* In other words, NW equals 2 to the w-th power */

#define TILE_WIDTH_ROW 2
#define TILE_WIDTH_COL 64
//#define TILE_WIDTH 2
#define TILE_DEPTH 2

#define BUFFER_SIZE 256

//#define DEBUG


//#define BLOCKSIZE 16
//#define BLOCKSIZEMINUS1 15
#define BLOCKSIZE 4
#define BLOCKSIZEMINUS1 3

//#define USELOOPUNROLLING 1  
#define AVOIDBANKCONFLICTS 0    //this just runs faster :X

__shared__ uint8_t gflog[512];
__shared__ uint8_t gfexp[512];

//__global__ int setup_tables(int w, uint8_t *gflog, uint8_t *gfexp)
//__device__ int setup_tables(int w, uint8_t *gflog, uint8_t *gfexp)
__device__ int setup_tables(int w)
{
	unsigned int b;
   	unsigned int log;
	unsigned int x_to_w;
	unsigned int prim_poly;
//	unsigned int r;
//	unsigned int x;
//	unsigned int y;

	unsigned int prim_poly_4 = 023;
	unsigned int prim_poly_8 = 0435;
	//uint8_t prim_poly_8 = 285;
	unsigned int prim_poly_16 = 0210013;
	switch(w) 
	{
		case 4: prim_poly = prim_poly_4; break;
		case 8: prim_poly = prim_poly_8; break;
		case 16: prim_poly = prim_poly_16; break;
		default: return -1;
	}
	x_to_w = 1 << w;
	b = 1;
//	r = 0;
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
		gfexp[log] = (uint8_t) b;
		b = b << 1;
		if (b & x_to_w) 
		{
			b = b ^ prim_poly;
		}
	}
	return 0;
}

__device__ uint8_t gf_add(uint8_t a, uint8_t b)
{
	return a^b;
}

__device__ uint8_t gf_sub(uint8_t a, uint8_t b)
{
	return gf_add(a, b);
}

__device__ uint8_t gf_mul(uint8_t a, uint8_t b)
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
	return gfexp[sum_log];
}

__device__ uint8_t gf_mul(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
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
	return gfexp[sum_log];
}

uint8_t gf_mul_bit(uint8_t a, uint8_t b)
{
	uint8_t sum_log;
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

__device__ uint8_t gf_mul_bit(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
{
	uint8_t sum_log;
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

__device__ uint8_t gf_div(uint8_t a, uint8_t b)
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
	return gfexp[diff_log];
}

__device__ uint8_t gf_div(uint8_t a, uint8_t b, uint8_t *gflog, uint8_t *gfexp)
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
	return gfexp[diff_log];
}

__device__ uint8_t gf_pow(uint8_t a, uint8_t power)
{
	int pow_log = (gflog[a] * power) % (NW-1);
	return gfexp[pow_log];
}

__device__ uint8_t gf_pow(uint8_t a, uint8_t power, uint8_t *gflog, uint8_t *gfexp)
{
	int pow_log = (gflog[a] * power) % (NW-1);
	return gfexp[pow_log];
}


// C=AB
// A: nxp
// B: pxm
// C: nxm
__device__ void matrix_mul(unsigned char *A, unsigned char *B, unsigned char *C, int n, int p, int m)
{
	__shared__ int rowVector[TILE_WIDTH_ROW][TILE_DEPTH];
	__shared__ int colVector[TILE_DEPTH][TILE_WIDTH_COL];
	__shared__ int product[TILE_WIDTH_ROW][TILE_WIDTH_COL];

	int bx = blockIdx.x;
   	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int row;
	int col;
	int px;
	int py;	

	setup_tables(8);
	__syncthreads();

	for(py=ty; py<TILE_WIDTH_ROW; py+=blockDim.y)
	{
		for(px=tx; px<TILE_WIDTH_COL; px+=blockDim.x)
		{
			row = by*TILE_WIDTH_ROW+py;
			col = bx*TILE_WIDTH_COL+px;
			product[py][px] = 0;
			__syncthreads();
		
			for(int i=0; i<(int)(ceil((float)p/TILE_DEPTH)); i++)
			{
				for(int j=tx; j<TILE_DEPTH; j+=blockDim.x)
				{
					rowVector[py][j] = A[row*p+i*TILE_DEPTH+j];
				}
				for(int j=ty; j<TILE_DEPTH; j+=blockDim.y)
				{		
					colVector[j][px] = B[col+(i*TILE_DEPTH+j)*m];
				}
				__syncthreads();
		
				for(int j=0; j<TILE_DEPTH; j++)
				{
					product[py][px] ^= gf_mul(rowVector[py][j], colVector[j][px]);
//					dist[py][px] = gf_add(dist[py][px], gf_mul(rowVector[py][j], colVector[j][px]));
				}
				__syncthreads();
			}
			C[row*m+col] = product[py][px];
		}
	}
	/*
	int i;
	int j;
	int k;
	setup_tables(8);
	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
		{
			for(k=0; k<p; k++)
			{
				C[i*m+j] = gf_add( C[i*m+j], gf_mul( A[i*p+k], B[k*m+j] ) );
			}
		}
	}
	*/
}

__global__ void decode_chunk(unsigned char *dataChunk, unsigned char *parityCoeff, unsigned char *codeChunk, int nativeBlockNum, int parityBlockNum, int chunkSize)
{
	matrix_mul(parityCoeff, codeChunk, dataChunk, nativeBlockNum, nativeBlockNum, chunkSize);
}

/******************************************************************************
                           AWAKE THE GPU CARD KERNEL

                                    works
*******************************************************************************/
__global__ void initKernel(){}


/******************************************************************************
                        ROUND UP - AIDING FUNCTION

*******************************************************************************/
int roundUp ( int n, int d )
{
	return n/d + (n%d != 0);
}

/******************************************************************************
                        MINIMUM - AIDING FUNCTION

*******************************************************************************/
int minimo(     int a,
                int b )
{
    return ( a < b )? a : b;

}

/******************************************************************************
                   SWITCH ROWS IF NECESSARY - KERNEL 

            Needs only a "row" of threads, block can be linear!

                                works                                    
*******************************************************************************/
__global__ void switchRows( uint8_t *matrix,
                            uint8_t *result,
                            int index, 
                            int rowToSwitch,
                            int lda )
{
//    int y = threadIdx.y + LINEAR_BLOCK_SIZE * blockIdx.y;
    int y = threadIdx.y + blockDim.y * blockIdx.y;
    uint8_t tmp_m, tmp_r;

    if ( y < lda )
    {
        tmp_m = matrix[ IDC2D( index, y, lda ) ];
        matrix[ IDC2D( index, y, lda ) ] = matrix[ IDC2D( rowToSwitch, y, lda ) ];
        matrix[ IDC2D( rowToSwitch, y, lda ) ] = tmp_m; 

        tmp_r = result[ IDC2D( index, y, lda ) ];
        result[ IDC2D( index, y, lda ) ] = result[ IDC2D( rowToSwitch, y, lda ) ];
        result[ IDC2D( rowToSwitch, y, lda ) ] = tmp_r; 
    }
} 
/******************************************************************************
                    NORMALIZE THE PIVOT ROW KERNEL

            Needs only a "column" of threads, block can be linear!

                                works
*******************************************************************************/
__global__ void normalizePivotRow(  uint8_t *matrix,
                                    uint8_t *result,
                                    int index,
                                    int lda )
{
// Variable declaration     
    int ty = threadIdx.y;                           // Position of each thread inside the block
//    int y = ty + LINEAR_BLOCK_SIZE * blockIdx.y;   // Position of each thread inside the matrix
    int y = ty + blockDim.y * blockIdx.y;   // Position of each thread inside the matrix

    __shared__ uint8_t pivotValue;        // Each "block" will have to load it!

	setup_tables(8);
	__syncthreads();

    if ( y < lda )      // CAN'T ACCESS TO VALUES OUTSIDE THE MATRIX! ONLY APPLIES FOR #threads > LINEAR_BLOCK_SIZE 
    {

    // First thread of each block loads the pivotValue
        if ( ty == 0 )
            pivotValue = matrix[ IDC2D( index, index, lda) ];
        __syncthreads();
    
    // Normalize the pivot row!
    // Every thread divides the element of its position with the pivotValue
        matrix[ IDC2D( index, y, lda )] = gf_div(matrix[ IDC2D( index, y, lda )], pivotValue);
        result[ IDC2D( index, y, lda )] = gf_div(result[ IDC2D( index, y, lda )], pivotValue);


    }

}




/******************************************************************************
         MODIFIED GAUSSIAN ELIMINATION (MGE) KERNEL - LINEAR VERSION

        As the name implies, its grid is formed out of linear blocks

*******************************************************************************/

__global__ void linearMge(  uint8_t *matrix,
                            uint8_t *result,
                            int index,
                            int lda )
{
    int ty = threadIdx.y;
    int x =  blockIdx.x;
//    int y = ty + blockIdx.y * LINEAR_BLOCK_SIZE;
    int y = ty + blockIdx.y * blockDim.y;

    __shared__ uint8_t multColumn[ LINEAR_BLOCK_SIZE ];

    __shared__ uint8_t matrixPivotValue;
    __shared__ uint8_t matrixRow[ LINEAR_BLOCK_SIZE ];
    
    __shared__ uint8_t resultPivotValue;
    __shared__ uint8_t resultRow[ LINEAR_BLOCK_SIZE];

    uint8_t newMatrixValue;
    uint8_t newResultValue;

	setup_tables(8);
	__syncthreads();

    if ( y < lda )
    {
    // Each block loads the value of the pivot Row to be substracted
        if ( ty == 0 )
        {
            matrixPivotValue = matrix[ IDC2D( index, x, lda )];
            resultPivotValue = result[ IDC2D( index, x, lda )];
        }
        multColumn[ ty ] = matrix[ IDC2D( y, index, lda )];
        
        matrixRow[ ty ] = matrix[ IDC2D( y, x, lda )]; 
        resultRow[ ty ] = result[ IDC2D( y, x, lda )]; 
        __syncthreads();

    // Now, the substraction, please :)
        if ( y!= index )
        {
            newMatrixValue = matrixRow[ty] ^ gf_mul(multColumn[ty], matrixPivotValue);
            newResultValue = resultRow[ty] ^ gf_mul(multColumn[ty], resultPivotValue);
        
        // Copy to the matrix
            matrix[ IDC2D( y, x, lda) ] = newMatrixValue;
            result[ IDC2D( y, x, lda) ] = newResultValue;
        }
    }
}



/******************************************************************************
                        CREATE RESULT KERNEL 

    Needs a complete grid of threads for ALL of the matrix positions.
    
    * SQUARE_BLOCK_SIZE and squareGrid used!

*******************************************************************************/

__global__ void createResult(   uint8_t *result,
                                int lda )
{

// Variable declaration 
    int x = blockIdx.x * SQUARE_BLOCK_SIZE + threadIdx.x;;  // Position identifier inside the matrix
    int y = blockIdx.y * SQUARE_BLOCK_SIZE + threadIdx.y;  // Position identifier inside the matrix


// Creation of the result matrix ( an identity matrix in first instance )
    if ( x == y )
        result[ IDC2D( y, x, lda) ] = 1;  // 1.0 in the diagonal
    else
        result[ IDC2D( y, x, lda) ] = 0;  // 0.0 elsewhere
}

/******************************************************************************
              GET THE MAXIMUM INDEX OF THE COLUMN - HOST CODE

                                    works
*******************************************************************************/
int getOkIndex( uint8_t *column,
                    int index,
                    int lda )
{
    int okIndex = -1;
    int i = index;
    uint8_t threshold = 0;

    while ( okIndex == -1 && i < lda )
    {
        okIndex = ( column[i] > threshold )? i : -1;        
        i++;
    }
    return okIndex;
}


/******************************************************************************
                MODIFIED GAUSSIAN ELIMINATION (MGE) - HOST CODE

*******************************************************************************/

void invert_matrix( uint8_t *matrix_dev, uint8_t *result_dev, int lda )
{
// Variable declaration
//    uint8_t *matrix_dev;
//    uint8_t *result_dev;

    uint8_t matrixColumn[lda];

    int sizeMat = lda * lda * sizeof( uint8_t );
    int sizeCol = lda * sizeof( uint8_t );

    int index;
    int okIndex;


// CUDA variable reservation
//    cudaMalloc( (void**) &matrix_dev, sizeMat );
//    cudaMemcpy( matrix_dev, input, sizeMat, cudaMemcpyHostToDevice );
//    cudaMalloc( (void**) &result_dev, sizeMat );

// setup for createResult kernel (cr)
    dim3 crGrid( roundUp( lda, SQUARE_BLOCK_SIZE ), roundUp( lda, SQUARE_BLOCK_SIZE ) );    // GRID
    dim3 crBlock( minimo( lda, SQUARE_BLOCK_SIZE ), minimo( lda, SQUARE_BLOCK_SIZE ) );     // BLOCK

// setup for normalizePivotRow kernel (npr) GRID
    dim3 nprGrid( 1, roundUp( lda, LINEAR_BLOCK_SIZE) );
 
// setup for linearMGE kernel (lmge) GRID
    dim3 lmgeGrid( lda, roundUp( lda, LINEAR_BLOCK_SIZE));

// setup for normalizePivotRow kernel (npr) AND linearMGE kernel BLOCK
    dim3 linearBlock( 1, minimo(lda, LINEAR_BLOCK_SIZE)); 

// Shows the setup if DISPLAY_SETTINGS macro is TRUE
    if ( DISPLAY_SETTINGS )
    {
        printf( "\nKernels Setup:\n" );

        printf( "\t> crGrid( %d, %d )\n", crGrid.x, crGrid.y );
        printf( "\t> crBlock( %d, %d )\n", crBlock.x, crBlock.y );    

        printf( "\t> nprGrid( %d, %d )\n", nprGrid.x, nprGrid.y );
        printf( "\t> nprBlock( %d, %d )\n", linearBlock.x, linearBlock.y );

        printf( "\t> lmgeGrid( %d, %d )\n", lmgeGrid.x, lmgeGrid.y );
        printf( "\t> lmgeBlock( %d, %d )\n\n", linearBlock.x, linearBlock.y );
    }

// Create my result, please :)
    createResult<<< crGrid, crBlock >>>( result_dev, lda );

// Launch my fucking algorithm, thx!

   for ( index = 0; index < lda; index ++ )
    {

    // Partial pivoting needed in bad cases...
        cudaMemcpy( matrixColumn, matrix_dev + index * lda , sizeCol , cudaMemcpyDeviceToHost );            

        okIndex = getOkIndex( matrixColumn, index, lda );

        if( DISPLAY )
             printf( "\t> okIndex = %d\n", okIndex );

        if ( okIndex != index )
            switchRows<<< nprGrid, linearBlock >>>( matrix_dev, result_dev, index, okIndex, lda );

    cudaThreadSynchronize();
    // Normalize the pivot row!
        normalizePivotRow<<< nprGrid, linearBlock>>>( matrix_dev, result_dev, index, lda );

    cudaThreadSynchronize();
    // MGE - linear version
        linearMge<<< lmgeGrid, linearBlock >>>( matrix_dev, result_dev, index, lda );
    }
    cudaThreadSynchronize();

//// Copy the damn result :)
//   cudaMemcpy( output, result_dev, sizeMat, cudaMemcpyDeviceToHost);
////    cudaMemcpy( output, matrix_dev, sizeMat, cudaMemcpyDeviceToHost);
//    cudaFree( matrix_dev ); cudaFree( result_dev );
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

void copy_matrix(uint8_t *src, uint8_t *des, int srcRowIndex, int desRowIndex, int rowSize)
{
	int i;
	for(i=0; i<rowSize; i++)
	{
		des[desRowIndex*rowSize+i] = src[srcRowIndex*rowSize+i];
	}
}

int main(int argc, char *argv[])
{
	int nativeBlockNum = 4;
	int parityBlockNum = 2;
	int chunkSize = 1;
	int totalSize;

	uint8_t *dataBuf;		//host
	uint8_t *codeBuf;		//host
	uint8_t *dataBuf_d;		//device
	uint8_t *codeBuf_d;		//device

	int dataSize;
	int codeSize;

	FILE *fp_in;
	FILE *fp_out;

	int totalMatrixSize;
	int matrixSize;
	uint8_t *totalEncodingMatrix;	//host
	uint8_t *encodingMatrix;	//host
	if( ( fp_in = fopen(".METADATA","rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}
	fscanf(fp_in, "%d", &totalSize);
	fscanf(fp_in, "%d %d", &parityBlockNum, &nativeBlockNum);
	chunkSize = (int) (ceil( (float)totalSize / nativeBlockNum )); 
	totalMatrixSize = nativeBlockNum * ( nativeBlockNum + parityBlockNum );
	totalEncodingMatrix = (uint8_t*) malloc( totalMatrixSize );
	matrixSize = nativeBlockNum * nativeBlockNum;
	encodingMatrix = (uint8_t*) malloc( matrixSize );
	for(int i =0; i<nativeBlockNum*(nativeBlockNum+parityBlockNum); i++)
	{
		fscanf(fp_in, "%d", totalEncodingMatrix+i);
	}

	dataSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	codeSize = nativeBlockNum*chunkSize*sizeof(uint8_t);
	dataBuf = (uint8_t*) malloc( dataSize );
	memset(dataBuf, 0, dataSize);
	codeBuf = (uint8_t*) malloc( codeSize);
	memset(codeBuf, 0, codeSize);
	cudaMalloc( (void **)&dataBuf_d, dataSize );
	cudaMemset(dataBuf_d, 0, dataSize);
	cudaMalloc( (void **)&codeBuf_d, codeSize );
	cudaMemset(codeBuf_d, 0, codeSize);

	if(argc == 2)
	{
		FILE *fp_conf;
		char input_file_name[20];
		int index;
		fp_conf = fopen(argv[1], "r");

		for(int i=0; i<nativeBlockNum; i++)
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
		for(int i=0; i<nativeBlockNum; i++)
		{
			char input_file_name[20];
			int index;
			printf("Please enter the file name of fragment:\n");
			scanf("%s", input_file_name);
			index = atoi(input_file_name+1);
			printf("#%dth fragment\n", index);

			copy_matrix(totalEncodingMatrix, encodingMatrix, index, i, nativeBlockNum);

			fp_in = fopen(input_file_name, "rb");
			fseek(fp_in, 0L, SEEK_SET);
			// this part can be process in parallel with computing inversed matrix
			fread(codeBuf+i*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
			fclose(fp_in);

		}
	}
	cudaMemcpy(codeBuf_d, codeBuf, codeSize, cudaMemcpyHostToDevice);

#ifdef DEBUG
	show_squre_matrix(encodingMatrix, nativeBlockNum);
#endif

/*
	if( ( fp_in = fopen("native_2","rb") ) == NULL )
	{
		printf("Can not open source file!\n");
		exit(0);
	}
	fseek(fp_in, 0L, SEEK_END);
	chunkSize = ftell(fp_in);
*/
#ifdef DEBUG
//	int matrixSize = nativeBlockNum*nativeBlockNum*sizeof(uint8_t);
//	uint8_t testMatrix[4][4] = {{0, 0, 1, 0}, {0, 1, 0, 0}, {1, 1, 1, 1}, {1, 2, 3, 4}};
//	uint8_t testMatrix[16] = {0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 3, 4};
//	uint8_t testMatrix[16] = {1, 1, 1, 1, 1, 2, 3, 4, 0, 0, 1, 0, 0, 0, 0, 1};
//	cudaMemcpy(encodingMatrix_d, testMatrix, matrixSize, cudaMemcpyHostToDevice);
#endif

	uint8_t *encodingMatrix_d;	//device
	uint8_t *decodingMatrix;	//host
	uint8_t *decodingMatrix_d;	//device
	decodingMatrix = (uint8_t*) malloc( matrixSize );
	cudaMalloc( (void **)&encodingMatrix_d, matrixSize );
	cudaMalloc( (void **)&decodingMatrix_d, matrixSize );
	cudaMemcpy(encodingMatrix_d, encodingMatrix, matrixSize, cudaMemcpyHostToDevice);

    invert_matrix( encodingMatrix_d, decodingMatrix_d, nativeBlockNum );

	cudaMemcpy(decodingMatrix, decodingMatrix_d, matrixSize, cudaMemcpyDeviceToHost);
	show_squre_matrix(decodingMatrix, nativeBlockNum);
/*
	fseek(fp_in, 0L, SEEK_SET);
	fread(codeBuf, sizeof(uint8_t), chunkSize, fp_in);
	fclose(fp_in);

	fp_in = fopen("native_1", "rb");
	fseek(fp_in, 0L, SEEK_SET);
	fread(codeBuf+chunkSize, sizeof(uint8_t), chunkSize, fp_in);
	fclose(fp_in);

	fp_in = fopen("code_0", "rb");
	fseek(fp_in, 0L, SEEK_SET);
	fread(codeBuf+2*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
	fclose(fp_in);

	fp_in = fopen("code_1", "rb");
	fseek(fp_in, 0L, SEEK_SET);
	fread(codeBuf+3*chunkSize, sizeof(uint8_t), chunkSize, fp_in);
	fclose(fp_in);
*/

//	int gridDimX = (int)(ceil((float)chunkSize/TILE_WIDTH));
//	int gridDimY = (int)(ceil((float)parityBlockNum/TILE_WIDTH));
//	dim3 block(TILE_WIDTH, TILE_WIDTH);
	int gridDimX = (int)( ceil((float)chunkSize / TILE_WIDTH_COL) );
	int gridDimY = (int)( ceil((float)nativeBlockNum / TILE_WIDTH_ROW) );
	dim3 grid(gridDimX, gridDimY);
	dim3 block(TILE_WIDTH_ROW, TILE_WIDTH_COL);
	decode_chunk<<<grid, block>>>(dataBuf_d, decodingMatrix_d, codeBuf_d, nativeBlockNum, parityBlockNum, chunkSize);
	cudaMemcpy(dataBuf, dataBuf_d, dataSize, cudaMemcpyDeviceToHost);

	char output_file_name[20];
	printf("Enter the name of the decoded file:\n");
	scanf("%s", output_file_name);
	fp_out = fopen(output_file_name, "wb");
	fwrite(dataBuf, sizeof(uint8_t), totalSize, fp_out);
	fclose(fp_out);

	cudaFree(decodingMatrix_d);
	cudaFree(dataBuf_d);
	cudaFree(codeBuf_d);

	free(decodingMatrix);
	free(dataBuf);
	free(codeBuf);

}

