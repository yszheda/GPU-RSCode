// =====================================================================================
// 
//       Filename:  GaloisFieldValue.cpp
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  2013/9/2 ä¸å 07:47:50
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
#include <cassert>
// using namespace std;

template < unsigned int gf_width >
constexpr unsigned int GaloisFieldValue<gf_width>::prim_poly_table[33];
// initialize static members
template < unsigned int gf_width >
bool GaloisFieldValue<gf_width>::full_table_ready = false;
template < unsigned int gf_width >
bool GaloisFieldValue<gf_width>::double_tables_ready = false;
template < unsigned int gf_width >
bool GaloisFieldValue<gf_width>::log_exp_tables_ready = false;
template < unsigned int gf_width >
bool GaloisFieldValue<gf_width>::log_exp_tables_v1_ready = false;
template < unsigned int gf_width >
bool GaloisFieldValue<gf_width>::log_exp_tables_v2_ready = false;
template < unsigned int gf_width >
bool GaloisFieldValue<gf_width>::log_exp_tables_v3_ready = false;

inline int index(int row, int col, int table_size)
{
		return row * table_size + col;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setupTables
// Description:  setup corresponding tables
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setupTables ()
{
		switch ( policy ) {
				case GaloisFieldPolicies::LOOP:	
						break;

				case GaloisFieldPolicies::FULLTABLE:	
						if (full_table_ready == false)
						{
								setupFullTables();
						}
						break;

				case GaloisFieldPolicies::DOUBLETABLES:	
						if (double_tables_ready == false)
						{
								setupDoubleTables();
						}
						break;

				case GaloisFieldPolicies::LOGEXPTABLES:
						if (log_exp_tables_ready == false)
						{
								setupLogExpTables();	
						}
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V1:
						if (log_exp_tables_v1_ready == false)
						{
								setupLogExpTablesV1();	
						}
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V2:
						if (log_exp_tables_v2_ready == false)
						{
								setupLogExpTablesV2();	
						}
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V3:
						if (log_exp_tables_v3_ready == false)
						{
								setupLogExpTablesV3();	
						}
						break;

				default:	
						break;
		}				/* -----  end switch  ----- */
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setupFullTables
// Description:  setup full tables
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setupFullTables ()
{
		static_assert(gf_width <= 30, "Cannot handle width which is larger than 30!");
		unsigned int prim_poly = prim_poly_table[gf_width];
		int field_size = 1 << gf_width;
		int table_size = field_size * field_size;
		gf_mul_table.reserve(table_size);
		gf_div_table.reserve(table_size);
		for (int i = 0; i < table_size; i++) {
				gf_mul_table.push_back(0);
				gf_div_table.push_back(0);
		}
		for (int i = 0; i < field_size; ++i)
		{
				for (int j = 0; j < field_size; ++j)
				{
						int prod = loopMul(i, j);
						gf_mul_table[index(i, j, table_size)] = prod;
						// safe for div table?
						// gf_div_table[index(prod, i, table_size)] = j;
						// gf_div_table[index(prod, j, table_size)] = i;
						if (j != 0)
						{
								gf_div_table[index(i, j, table_size)] = loopDiv(i, j);
						}
				}
		}
		full_table_ready = true;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setupDoubleTables
// Description:  setup double tables
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setupDoubleTables ()
{
		static_assert(gf_width <= 30, "Cannot handle width which is larger than 30!");
		unsigned int prim_poly = prim_poly_table[gf_width];
		int field_size = 1 << gf_width;
		int table_size = (1 << (gf_width / 2)) * field_size;
		gf_left_mul_table.reserve(table_size);
		gf_right_mul_table.reserve(table_size);
		for (int i = 0; i < table_size; i++) {
				gf_left_mul_table.push_back(0);
				gf_right_mul_table.push_back(0);
		}
		for (int i = 0; i < (field_size >> gf_width/2); ++i)
		{
				for (int j = 0; j < field_size; ++j)
				{
						gf_left_mul_table[index(i, j, table_size)] = loopMul(i << gf_width/2, j);
						gf_right_mul_table[index(i, j, table_size)] = loopMul(i, j);
						gf_left_div_table[index(i, j, table_size)] = loopDiv(i << gf_width/2, j);
						gf_right_div_table[index(i, j, table_size)] = loopDiv(i, j);
				}
		}
		double_tables_ready = true;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setupLogExpTables
// Description:  setup log and exp tables
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setupLogExpTables ()
{
		static_assert(gf_width <= 30, "Cannot handle width which is larger than 30!");
		unsigned int prim_poly = prim_poly_table[gf_width];
		int table_size = 1 << gf_width;
		int field_size = 1 << gf_width;
		gf_log_table.reserve(table_size);
		gf_exp_table.reserve(table_size);
		for (int i = 0; i < table_size; i++) {
				gf_log_table.push_back(0);
				gf_exp_table.push_back(0);
		}
		int exp = 1;
		for (int i = 0; i < field_size - 1; i++) {
				if (exp > field_size) {
						// exception?
						break;
				}
				gf_log_table[exp] = i;
				gf_exp_table[i] = exp;
				exp = exp << 1;
				if (exp & field_size) {
						exp = exp ^ prim_poly;
				}
		}
		log_exp_tables_ready = true;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setupLogExpTablesV1
// Description:  setup log and exp tables with improvement I
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setupLogExpTablesV1 ()
{
		static_assert(gf_width <= 30, "Cannot handle width which is larger than 30!");
		unsigned int prim_poly = prim_poly_table[gf_width];
		int field_size = 1 << gf_width;
		int gf_max_value = field_size - 1;
		int log_table_size = field_size;
		int exp_table_size = field_size + 1;
		gf_log_table.reserve(log_table_size);
		gf_exp_table.reserve(exp_table_size);
		for (int i = 0; i < log_table_size; i++) {
				gf_log_table.push_back(0);
		}
		for (int i = 0; i < exp_table_size; i++) {
				gf_log_table.push_back(0);
		}
		int exp = 1;
		for (int i = 0; i < field_size - 1; i++) {
				if (exp > field_size) {
						// exception?
						break;
				}
				gf_log_table[exp] = i;
				gf_exp_table[i] = exp;
				exp = exp << 1;
				if (exp & field_size) {
						exp = exp ^ prim_poly;
				}
		}
		gf_exp_table[gf_max_value] = gf_exp_table[0];
		log_exp_tables_v1_ready = true;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setupLogExpTablesV2
// Description:  setup log and exp tables with improvement II
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setupLogExpTablesV2 ()
{
		static_assert(gf_width <= 30, "Cannot handle width which is larger than 30!");
		unsigned int prim_poly = prim_poly_table[gf_width];
		int field_size = 1 << gf_width;
		int gf_max_value = field_size - 1;
		int log_table_size = field_size;
		int exp_table_size = gf_max_value * 2 - 2;
		gf_log_table.reserve(log_table_size);
		gf_exp_table.reserve(exp_table_size);
		for (int i = 0; i < log_table_size; i++) {
				gf_log_table.push_back(0);
		}
		for (int i = 0; i < exp_table_size; i++) {
				gf_log_table.push_back(0);
		}
		int exp = 1;
		for (int i = 0; i < field_size - 1; i++) {
				if (exp > field_size) {
						// exception?
						break;
				}
				gf_log_table[exp] = i;
				gf_exp_table[i] = exp;
				if (i < gf_max_value - 1)
				{
					gf_exp_table[i + gf_max_value] = exp;
				}
				exp = exp << 1;
				if (exp & field_size) {
						exp = exp ^ prim_poly;
				}
		}
		log_exp_tables_v2_ready = true;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  setupLogExpTablesV3
// Description:  setup log and exp tables with improvement III
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::setupLogExpTablesV3 ()
{
		static_assert(gf_width <= 30, "Cannot handle width which is larger than 30!");
		unsigned int prim_poly = prim_poly_table[gf_width];
		int field_size = 1 << gf_width;
		int gf_max_value = field_size - 1;
		int log_table_size = gf_max_value + 1;
		int exp_table_size = gf_max_value * 4 + 1;
		gf_log_table.reserve(log_table_size);
		gf_exp_table.reserve(exp_table_size);
		for (int i = 0; i < log_table_size; i++) {
				gf_log_table.push_back(0);
		}
		for (int i = 0; i < exp_table_size; i++) {
				gf_log_table.push_back(0);
		}
		gf_log_table[0] = 2 * gf_max_value;
		int exp = 1;
		for (int i = 0; i < field_size - 1; i++) {
				if (exp > field_size) {
						// exception?
						break;
				}
				gf_log_table[exp] = i;
				gf_exp_table[i] = exp;
				if (i < field_size - 1)
				{
					gf_exp_table[i + gf_max_value] = exp;
				}
				exp = exp << 1;
				if (exp & field_size) {
						exp = exp ^ prim_poly;
				}
		}
		log_exp_tables_v3_ready = true;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  destroyTables
// Description:  destroy all tables
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void 
GaloisFieldValue<gf_width>::destroyTables ()
{
		// completely clear vectors
		std::vector<int>().swap(gf_mul_table);
		std::vector<int>().swap(gf_div_table);
		std::vector<int>().swap(gf_left_mul_table);
		std::vector<int>().swap(gf_right_mul_table);
		std::vector<int>().swap(gf_left_div_table);
		std::vector<int>().swap(gf_right_div_table);
		std::vector<int>().swap(gf_log_table);
		std::vector<int>().swap(gf_exp_table);

		// reset flags
		full_table_ready = false;
		double_tables_ready = false;
		log_exp_tables_ready = false;
		log_exp_tables_v1_ready = false;
		log_exp_tables_v2_ready = false;
		log_exp_tables_v3_ready = false;
}




//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  default constructor
// Description:  default constructor
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>::GaloisFieldValue ( ): gf_value(0), policy(GaloisFieldPolicies::LOOP)
{
		this->setupTables();
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  constructor
// Description:  constructor
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>::GaloisFieldValue ( const int gf_value ): gf_value(gf_value)
{
		this->setupTables();
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  copy constructor
// Description:  copy constructor
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>::GaloisFieldValue ( const GaloisFieldValue<gf_width> &other ): gf_value(other.gf_value), policy(other.policy)
{
		this->destroyTables();
		this->setupTables();
}

//--------------------------------------------------------------------------------------
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
				policy = other.policy;

		}
		this->destroyTables();
		this->setupTables();
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
//      Method:  loopMul
// Description:  multiplication by loop-based method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
int
GaloisFieldValue<gf_width>::loopMul ( const int x, const int y )
{
		const int gf_max_value = (1 << gf_width) - 1;
		int result;
		while(y) {
				if (y & 1)
				{
						result ^= x;
				}
				x <<= 1;
				if (x & (1 << (gf_width-1)))
				{
						x = (x ^ prim_poly_table[gf_width]) & gf_max_value;
				}
				y >>= 1;
		}
		return result;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  loopDiv
// Description:  division by loop-based method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
int
GaloisFieldValue<gf_width>::loopDiv( const int x, const int y )
{
		if (x == 0)
		{
				return 0;
		}
		if (y == 0)
		{
				// exception?
				return -1;
		}
		int inverse = loopInv(y);
		return loopMul(x, inverse);
}


//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  loopInv
// Description:  inversion by loop-based method
// adapted from Jeasure
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
int
GaloisFieldValue<gf_width>::loopInv( const int x )
{
		if (x == 0)
		{
				// exception?
				return -1;
		}
		int mat[32];
		int inv[32];
		const int gf_max_value = (1 << gf_width) - 1;
		for (int i = 0; i < gf_width; ++i)
		{
				mat[i] = x;
				x <<= 1;
				if (x & (1 << (gf_width-1)))
				{
						x = (x ^ prim_poly_table[gf_width]) & gf_max_value;
				}
		}
		invertBinaryMatrix(mat, inv);
		return inv[0];
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  invertBinaryMatrix
// Description:  
// adapted from Jeasure
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
void
GaloisFieldValue<gf_width>::invertBinaryMatrix( int* mat, int* inv )
{
		int size = gf_width;
		for (int i = 0; i < size; ++i)
		{
				inv[i] = 1 << i;
		}
		for (int i = 0; i < size; ++i)
		{
				int j;
				if ((mat[i] & (1 << i)) == 0)
				{
						for (j = i+1; j < size && (mat[j] & (1 << i)) == 0; ++j) ;
						if (j == size)
						{
								// invertible
								exit(1);
						}
						std::swap(mat[i], mat[j]);
						std::swap(inv[i], inv[j]);
						
				}
				for (j = i+1; j < size; ++j)
				{
						if ((mat[j] & (1 << i)) != 0)
						{
								mat[j] ^= mat[i];
								inv[j] ^= inv[i];
						}
				}
		}
		for (int i = size - 1; i >= 0; --i)
		{
				for (int j = 0; j < i; ++j)
				{
						if (mat[j] & (1 << i))
						{
								inv[j] ^= inv[i];
						}
				}
		}
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  loopMulAssign
// Description:  mul assignment operator by loop-based method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::loopMulAssign ( const GaloisFieldValue<gf_width> &other )
{
		this->gf_value = loopMul(this->gf_value, other.gf_value);
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  fullTableMulAssign
// Description:  mul assignment operator by full-table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::fullTableMulAssign ( const GaloisFieldValue<gf_width> &other )
{
//		assert(full_table_ready == true);
		int field_size = 1 << gf_width;
		int table_size = field_size * field_size;
		this->gf_value = gf_mul_table[index(this->gf_value, other.gf_value, table_size)];
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  doubleTablesMulAssign
// Description:  mul assignment operator by double-table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::doubleTablesMulAssign ( const GaloisFieldValue<gf_width> &other )
{
//		assert(double_tables_ready == true);
		int mask = (1 << (gf_width/2)) - 1;
		int left_value = (this->gf_value >> (gf_width/2)) & mask;
		int right_value = (this->gf_value) & mask;
		int field_size = 1 << gf_width;
		int table_size = (1 << (gf_width / 2)) * field_size;
		this->gf_value = gf_left_mul_table[index(left_value, other.gf_value, table_size)] 
				^ gf_right_mul_table[index(right_value, other.gf_value, table_size)];
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesMulAssign
// Description:  mul assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesMulAssign ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_ready == true);
		if (this->gf_value == 0 || other.gf_value == 0) {
				gf_value = 0;
		} else {
			int gf_max_value = (1 << gf_width) - 1;
			int sum_log = gf_log_table[gf_value] + gf_log_table[other.gf_value];
			if (sum_log >= gf_max_value) {
					sum_log -= gf_max_value;
			}
			gf_value = gf_exp_table[sum_log];
		}
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesMulAssignV0
// Description:  mul assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesMulAssignV0 ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_ready == true);
		if (this->gf_value == 0 || other.gf_value == 0) {
				gf_value = 0;
		} else {
			int gf_max_value = (1 << gf_width) - 1;
			int sum_log = gf_log_table[gf_value] + gf_log_table[other.gf_value];
			gf_value = gf_exp_table[sum_log % gf_max_value];
		}
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesMulAssignV1
// Description:  mul assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesMulAssignV1 ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_v1_ready == true);
		if (this->gf_value == 0 || other.gf_value == 0) {
				gf_value = 0;
		} else {
			int gf_max_value = (1 << gf_width) - 1;
			int sum_log = gf_log_table[gf_value] + gf_log_table[other.gf_value];
			gf_value = gf_exp_table[sum_log & gf_max_value + sum_log >> gf_width];
		}
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesMulAssignV2
// Description:  mul assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesMulAssignV2 ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_v2_ready == true);
		if (this->gf_value == 0 || other.gf_value == 0) {
				gf_value = 0;
		} else {
			int sum_log = gf_log_table[gf_value] + gf_log_table[other.gf_value];
			gf_value = gf_exp_table[sum_log];
		}
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesMulAssignV3
// Description:  mul assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesMulAssignV3 ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_v3_ready == true);
		int sum_log = gf_log_table[gf_value] + gf_log_table[other.gf_value];
		gf_value = gf_exp_table[sum_log];
		return *this;
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
		switch ( policy ) {
				case GaloisFieldPolicies::LOOP:	
						loopMulAssign(other);
						break;

				case GaloisFieldPolicies::FULLTABLE:	
						fullTableMulAssign(other);
						break;

				case GaloisFieldPolicies::DOUBLETABLES:	
						doubleTablesMulAssign(other);
						break;

				case GaloisFieldPolicies::LOGEXPTABLES:
						logExpTablesMulAssign(other);	
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V1:
						logExpTablesMulAssignV1(other);
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V2:
						logExpTablesMulAssignV2(other);
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V3:
						logExpTablesMulAssignV3(other);
						break;

				default:	
						return *this;
						break;
		}				/* -----  end switch  ----- */
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  loopDivAssign
// Description:  div assignment operator by loop-based method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::loopDivAssign ( const GaloisFieldValue<gf_width> &other )
{
		this->gf_value = loopDiv(this->gf_value, other.gf_value);
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  fullTableDivAssign
// Description:  div assignment operator by full-table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::fullTableDivAssign ( const GaloisFieldValue<gf_width> &other )
{
//		assert(full_table_ready == true);
		int field_size = 1 << gf_width;
		int table_size = field_size * field_size;
		this->gf_value = gf_div_table[index(this->gf_value, other.gf_value, table_size)];
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  doubleTablesDivAssign
// Description:  div assignment operator by double-table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::doubleTablesDivAssign ( const GaloisFieldValue<gf_width> &other )
{
//		assert(double_tables_ready == true);
		int mask = (1 << (gf_width/2)) - 1;
		int left_value = (this->gf_value >> (gf_width/2)) & mask;
		int right_value = (this->gf_value) & mask;
		int field_size = 1 << gf_width;
		int table_size = (1 << (gf_width / 2)) * field_size;
		this->gf_value = gf_left_div_table[index(left_value, other.gf_value, table_size)] 
				^ gf_right_div_table[index(right_value, other.gf_value, table_size)];
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesDivAssign
// Description:  div assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesDivAssign ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_ready == true);
		static_assert(other.gf_value != 0, "Cannot divide by zero!");
		if (this->gf_value == 0) {
				gf_value = 0;
		} else {
			int gf_max_value = (1 << gf_width) - 1;
			int diff_log = gf_log_table[gf_value] - gf_log_table[other.gf_value];
			if (diff_log < 0) {
					diff_log += gf_max_value;
			}
			gf_value = gf_exp_table[diff_log];
		}
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesDivAssignV1
// Description:  div assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesDivAssignV1 ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_v1_ready == true);
		logExpTablesDivAssign(other);
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesDivAssignV2
// Description:  div assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesDivAssignV2 ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_v2_ready == true);
		logExpTablesDivAssign(other);
		return *this;
}

//--------------------------------------------------------------------------------------
//       Class:  GaloisFieldValue
//      Method:  logExpTablesDivAssignV3
// Description:  div assignment operator by log&exp table method
//--------------------------------------------------------------------------------------
template < unsigned int gf_width >
GaloisFieldValue<gf_width>&
GaloisFieldValue<gf_width>::logExpTablesDivAssignV3 ( const GaloisFieldValue<gf_width> &other )
{
//		assert(log_exp_tables_v3_ready == true);
		const int gf_max_value = (1 << gf_width) - 1;
		int diff_log = gf_log_table[gf_value] + gf_max_value - gf_log_table[other.gf_value];
		gf_value = gf_exp_table[diff_log];
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
		switch ( policy ) {
				case GaloisFieldPolicies::LOOP:	
						loopDivAssign(other);
						break;

				case GaloisFieldPolicies::FULLTABLE:	
						fullTableDivAssign(other);
						break;

				case GaloisFieldPolicies::DOUBLETABLES:	
						doubleTablesDivAssign(other);
						break;

				case GaloisFieldPolicies::LOGEXPTABLES:
						logExpTablesDivAssign(other);	
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V1:
						logExpTablesDivAssignV1(other);
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V2:
						logExpTablesDivAssignV2(other);
						break;

				case GaloisFieldPolicies::LOGEXPTABLES_V3:
						logExpTablesDivAssignV3(other);
						break;

				default:	
						return *this;
						break;
		}				/* -----  end switch  ----- */
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
		int gf_max_value = (1 << gf_width) - 1;
		int pow_log = (gf_log_table[gf_value] * other) % gf_max_value;
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

