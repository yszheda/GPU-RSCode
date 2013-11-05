/*
 * =====================================================================================
 *
 *       Filename:  ReedSolomon.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/05/13 17:35:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuai YUAN(yszheda AT gmail.com), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "Matrix.hpp"
#include "GaloisFieldValue.hpp"

template < unsigned int gf_width = 8 >
void encodeBlocks ( const int native_block_num, 
		const int parity_block_num, 
		const int block_size, 
		const GaloisFieldValue<gf_width>& data_blocks, 
		const GaloisFieldValue<gf_width>& parity_coeff, 
		GaloisFieldValue<gf_width>& code_blocks )
{
	matrixMul< GaloisFieldValue<gf_width> >(parity_block_num, native_block_num, block_size, parity_coeff, data_blocks, code_blocks);
}		/* -----  end of template function encodeBlocks  ----- */

template < unsigned int gf_width = 8 >
void decodeBlocks ( const int native_block_num, 
		const int block_size, 
		GaloisFieldValue<gf_width>& data_blocks, 
		const GaloisFieldValue<gf_width>& parity_coeff, 
		const GaloisFieldValue<gf_width>& code_blocks )
{
	matrixMul< GaloisFieldValue<gf_width> >(native_block_num, native_block_num, block_size, parity_coeff, code_blocks, data_blocks);
}		/* -----  end of template function decodeBlocks  ----- */

