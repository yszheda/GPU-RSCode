// =====================================================================================
// 
//       Filename:  GaloisFieldValue.cpp
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  2013/9/2 下午 07:47:50
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Shuai YUAN (galoisplusplus), yszheda AT gmail DOT com
//        Company:  
// 
// =====================================================================================
#include "GaloisFieldValue.hpp"
#include <iostream>
#include <vector>
// using namespace std;

template < unsigned int gf_width >
constexpr unsigned int GaloisFieldValue<gf_width>::prim_poly_table[33];
//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setup_tables
// Description:  setup log and exp tables
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setup_tables ()
{
		switch ( policy ) {
				case POLYMULT:	
						break;

				case FULLTABLE:	
						setup_full_tables();
						break;

				case LOGEXPTABLES:
						setup_log_exp_tables();	
						break;

				default:	
						break;
		}				/* -----  end switch  ----- */
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setup_tables
// Description:  setup log and exp tables
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setup_tables ()
{
		static_assert(gf_width <= 30, "Cannot handle width which is larger than 30!");
		unsigned int prim_poly = prim_poly_table[gf_width];
		int table_size = (1 << gf_width) - 1;
		int exp = 1;
		for (unsigned int i = 0; i < table_size; i++) {
				if (exp > table_size + 1) {
						break;
				}
				gf_log_table[exp] = i;
				gf_exp_table[i] = exp;
				exp = exp << 1;
				if (exp & table_size + 1) {
						exp = exp ^ prim_poly;
				}
		}
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  operator =
// Description:  assignment operator
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::operator = ( const GaloisFieldValue &other )
{
		if ( this != &other ) {
				gf_value = other.gf_value;
		}
		return *this;
}  // -----  end of method GaloisFieldValue::operator =  (assignment operator)  -----

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  operator +=
// Description:  add assignment operator
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::operator += ( const GaloisFieldValue<gf_width> &other )
{
		gf_value ^= other.gf_value; 	// XOR
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  operator -=
// Description:  sub assignment operator
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::operator -= ( const GaloisFieldValue<gf_width> &other )
{
//		operator+=(other);
		return (*this += other);
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  operator *=
// Description:  mul assignment operator
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::operator *= ( const GaloisFieldValue<gf_width> &other )
{
		if (this.gf_value == 0 || other.gf_value == 0) {
				gf_value = 0;
		} else {
			int gf_max_num = (1 << gf_width) - 1;
			int sum_log = gf_log_table[gf_value] + gf_log_table[other.gf_value];
			if (sum_log >= gf_max_num) {
					sum_log -= gf_max_num;
			}
			gf_value = gf_exp_table[sum_log];
		}
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  operator /=
// Description:  div assignment operator
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::operator /= ( const GaloisFieldValue<gf_width> &other )
{
		static_assert(other.gf_value != 0, "Cannot divide by zero!");
		if (this.gf_value == 0) {
				gf_value = 0;
		} else {
			int gf_max_num = (1 << gf_width) - 1;
			int diff_log = gf_log_table[gf_value] - gf_log_table[other.gf_value];
			if (diff_log < 0) {
					diff_log += gf_max_num;
			}
			gf_value = gf_exp_table[diff_log];
		}
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  operator ^=
// Description:  pow assignment operator
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::operator ^= ( const int &other )
{
		int gf_max_num = (1 << gf_width) - 1;
		int pow_log = (gf_log_table[gf_value] * other) % gf_max_num;
		gf_value = gf_exp_table[pow_log];
		return *this;
}

// ===  FUNCTION  ======================================================================
//         Name:  operator <<
//  Description:  
// =====================================================================================
template < unsigned int gf_width >
std::ostream& operator << ( std::ostream &os, const GaloisFieldValue<gf_width> &value )
{
		os << value.gf_value;
		return os;
}

// ===  FUNCTION  ======================================================================
//         Name:  operator +
//  Description:  add operator 
// =====================================================================================
template < unsigned int gf_width >
GaloisFieldValue<gf_width> operator + ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs )
{
		GaloisFieldValue<gf_width> result = lhs;
		result += rhs;
		return result;
}

// ===  FUNCTION  ======================================================================
//         Name:  operator -
//  Description:  sub operator 
// =====================================================================================
template < unsigned int gf_width >
GaloisFieldValue<gf_width> operator - ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs )
{
		GaloisFieldValue<gf_width> result = lhs;
		result -= rhs;
		return result;
}

// ===  FUNCTION  ======================================================================
//         Name:  operator *
//  Description:  mul operator 
// =====================================================================================
template < unsigned int gf_width >
GaloisFieldValue<gf_width> operator * ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs )
{
		GaloisFieldValue<gf_width> result = lhs;
		result *= rhs;
		return result;
}

// ===  FUNCTION  ======================================================================
//         Name:  operator /
//  Description:  div operator 
// =====================================================================================
template < unsigned int gf_width >
GaloisFieldValue<gf_width> operator / ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs )
{
		GaloisFieldValue<gf_width> result = lhs;
		result /= rhs;
		return result;
}

// ===  FUNCTION  ======================================================================
//         Name:  operator ^
//  Description:  add operator 
// =====================================================================================
template < unsigned int gf_width >
GaloisFieldValue<gf_width> operator ^ ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs )
{
		GaloisFieldValue<gf_width> result = lhs;
		result ^= rhs;
		return result;
}

