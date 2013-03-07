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

typedef unsigned char gf;

#define W 4
#define NW (1 << W) /* In other words, NW equals 2 to the w-th power */


gf gflog[16];
gf gfexp[16];
gf gfmul_16[16][16];

// int setup_tables(int w)
void setup_tables()
{
	unsigned int b;
   	unsigned int log;
	unsigned int x_to_w;
	unsigned int prim_poly;

	unsigned int prim_poly_4 = 023;
	unsigned int prim_poly_8 = 0435;
	unsigned int prim_poly_16 = 0210013;
	/* 
	switch(w) 
	{
		case 4: prim_poly = prim_poly_4; break;
		case 8: prim_poly = prim_poly_8; break;
		case 16: prim_poly = prim_poly_16; break;
		default: return -1;
	}
	x_to_w = 1 << w;
	*/

	prim_poly = prim_poly_4;
	x_to_w = 1 << 4;

	b = 1;
	for (log = 0; log < x_to_w-1; log++) 
	{
		if(b > x_to_w) break;
		gflog[b] = (gf) log;
		gfexp[log] = (gf) b;
		b = b << 1;
		if (b & x_to_w) 
		{
			b = b ^ prim_poly;
		}
	}
//	return 0;
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
	sum_log = gflog[a] + gflog[b];
	if (sum_log >= NW-1)
	{	
		sum_log -= NW-1;
	}
	return gfexp[sum_log];
}

gf gf_div(gf a, gf b)
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
	diff_log = gflog[a] - gflog[b];
	if (diff_log < 0)
	{	
		diff_log += NW-1;
	}
	return gfexp[diff_log];
}

gf gf_pow(gf a, gf power)
{
	int pow_log = (gflog[a] * power) % (NW-1);
	return gfexp[pow_log];
}

void gen_mult_table()
{
	int i;
	int j;
	for(i=0; i<16; i++)
	{
		for(j=0; j<16; j++)
		{
			gfmul_16[i][j] = gf_mul(i, j);
		}
	}
}

void show_mult_table()
{
	int i;
	int j;
	for(i=0; i<16; i++)
	{
		for(j=0; j<16; j++)
		{
			printf("%d ", gfmul_16[i][j]);
		}
		printf("\n");
	}
}

void print_exp_table(FILE *fp)
{
	int i;
	fprintf(fp, "__shared__ uint8_t gfexp[16] = {");
	for(i=0; i<15; i++)
	{
		fprintf(fp, " %d, ", gfexp[i]);
	}
	fprintf(fp, " %d }; \n", gfexp[15]);
}

void print_log_table(FILE *fp)
{
	int i;
	fprintf(fp, "__shared__ uint8_t gflog[16] = {");
	for(i=0; i<15; i++)
	{
		fprintf(fp, " %d, ", gflog[i]);
	}
	fprintf(fp, " %d }; \n", gflog[15]);
}



void print_mult_table(FILE *fp)
{
	int i;
	int j;
	fprintf(fp, "__shared__ uint8_t gfmul_16[16][16] = {");
	for(i=0; i<16; i++)
	{
		fprintf(fp, " {");
		for(j=0; j<16; j++)
		{
			fprintf(fp, " %d ", gfmul_16[i][j]);
			if(j == 15)
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


int main()
{
	setup_tables();
	gen_mult_table();
	show_mult_table();

	FILE *fp;
	fp = fopen("gf16.h", "w");
	fprintf(fp, "#ifdef _GF16_H_ \n");
	fprintf(fp, "#define _GF16_H_ \n");
	print_exp_table(fp);
	print_log_table(fp);
	print_mult_table(fp);
	fprintf(fp, "#endif");

}

