/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/25/2012 04:32:42 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuai YUAN (yszheda AT gmail.com), 
 *        Company:  
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "encode.h"
#include "decode.h"

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
/*	
	if(op == encode)
	{
		file = argv[4];
		encode_file(file, nativeBlockNum, parityBlockNum);
	}

	if(op == decode)
	{
		file = argv[4];
		decode_file(file, nativeBlockNum, parityBlockNum);
	}
*/

	switch(op)
	{
		case encode:
			inFile = argv[4];
			encode_file(inFile, nativeBlockNum, parityBlockNum);
			break;

		case decode:
			/*
			if(argc > 4)
			{
				file = argv[4];
				decode_file(file, nativeBlockNum, parityBlockNum);
			}
			else
			{
				decode_file(NULL, nativeBlockNum, parityBlockNum);
			}
			*/
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
