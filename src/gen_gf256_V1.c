/*
 * =====================================================================================
 *
 *       Filename:  gen_gf16.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/06/2013 08:35:49 PM
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
#include <math.h>

typedef uint8_t gf;

const int gf_width = 8;
const int field_size = 1 << 8;

gf gflog[256];
gf gfexp[256];

int setup_tables()
{
//	const int gf_width = 8;
//	const int field_size = 1 << gf_width;
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
	gfexp[field_size - 1] = gfexp[0];
	return 0;
}

gf gf_add(gf a, gf b)
{
	return a^b;
}

gf gf_sub(gf a, gf b)
{
	return gf_add(a, b);
}


gf gf_mul(gf a, gf b)
{
	int sum_log;
	if (a == 0 || b == 0)
	{
		return 0;
	}
	sum_log = (gflog[a] + gflog[b]) % (field_size - 1);
	return gfexp[sum_log];
}

gf gf_div(gf a, gf b)
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

gf gf_pow(gf a, int power)
{
	int pow_log = (gflog[a] * power) % (field_size - 1);
	return gfexp[pow_log];
}

/*  
void gen_mult_table()
{
	int i;
	int j;
	for(i = 0; i < field_size; i++)
	{
		for(j = 0; j < field_size; j++)
		{
			gfmul[i][j] = gf_mul(i, j);
		}
	}
}

void show_mult_table()
{
	int i;
	int j;
	for(i = 0; i < field_size; i++)
	{
		for(j = 0; j < field_size; j++)
		{
			printf("%d ", gfmul[i][j]);
		}
		printf("\n");
	}
}
*/

void print_exp_table(FILE *fp)
{
	int i;
	fprintf(fp, "__const__ uint8_t gfexp[%d] = {", field_size);
	for(i = 0; i < field_size - 1; i++)
	{
		fprintf(fp, " %d, ", gfexp[i]);
	}
	fprintf(fp, " %d }; \n", gfexp[field_size - 1]);
}

void print_log_table(FILE *fp)
{
	int i;
	fprintf(fp, "__const__ uint8_t gflog[%d] = {", field_size);
	for(i = 0; i < field_size - 1; i++)
	{
		fprintf(fp, " %d, ", gflog[i]);
	}
	fprintf(fp, " %d }; \n", gflog[field_size - 1]);
}

/*
void print_mult_table(FILE *fp)
{
	int i;
	int j;
	fprintf(fp, "__const__ uint8_t gfmul[%d][%d] = {", field_size, field_size);
	for(i = 0; i < field_size; i++)
	{
		fprintf(fp, " {");
		for(j = 0; j < field_size; j++)
		{
			fprintf(fp, " %d ", gfmul[i][j]);
			if(j == field_size - 1)
			{
				fprintf(fp, "}, \\ \n");
			}
			else
			{
				fprintf(fp, ", ");
			}
		}
//		fprintf(fp, "},\n");
	}
	fprintf(fp, "}; \n");
}
*/

int main()
{
	setup_tables();
//	gen_mult_table();
//	show_mult_table();

	FILE *fp;
	fp = fopen("gf256_V1.h", "w");
	fprintf(fp, "#ifdef _GF256_H_ \n");
	fprintf(fp, "#define _GF256_H_ \n");
	print_exp_table(fp);
	print_log_table(fp);
//	print_mult_table(fp);
	fprintf(fp, "#endif");
	fclose(fp);

	return 0;
}

