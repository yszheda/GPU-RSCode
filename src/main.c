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
#include <getopt.h>
#include <assert.h>
#include "encode.h"
#include "decode.h"

void show_help_info()
{
	printf("Usage:\n");
	printf("[-h]: show usage information\n");
	printf("Encode: [-k|-K nativeBlockNum] [-n|-N totalBlockNum] [-e|-E fileName]\n");
	printf("Decode: [-d|-D] [-k|-K nativeBlockNum] [-n|-N totalBlockNum] \n\t [-i|-I originalFileName] [-c|-C config] [-o|-O output]\n");
	printf("For encoding, the -k, -n, and -e options are all necessary.\n");
	printf("For decoding, the -d, -i, and -c options are all necessary.\n");
	printf("If the -o option is not set, the original file name will be chosen as the output file name by default.\n");
	printf("Performance-tuning Options:\n");
	printf("[-p]: set maxmimum blockDimX\n");
	printf("[-s]: set stream number\n");
	exit(0);
}

int main(int argc, char *argv[])
{
	int nativeBlockNum = 0;
	int parityBlockNum = 0;
	int totalBlockNum = 0;
	int gridDimXSize = 0;
	int streamNum = 1;
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
	while((option = getopt(argc, argv, "Ss:Pp:Kk:Nn:Ee:Ii:Cc:Oo:Ddh")) != -1) {
		switch (option) {
			case 'S':	
			case 's':	
				streamNum = (int) atoi(optarg);
				break;

			case 'P':	
			case 'p':	
				gridDimXSize = (int) atoi(optarg);
				break;

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
	
	switch ( op ) {
		case encode:	
			assert(nativeBlockNum != 0);
			assert(totalBlockNum != 0);
			parityBlockNum = totalBlockNum - nativeBlockNum;
			encode_file(inFile, nativeBlockNum, parityBlockNum, gridDimXSize, streamNum);
			break;

		case decode:	
			assert(inFile != NULL);
			assert(confFile != NULL);
			assert(outFile != NULL);
			decode_file(inFile, confFile, outFile, gridDimXSize, streamNum);
			break;

		default:	
			break;
	}		/* -----  end switch  ----- */

	return 0;

}
