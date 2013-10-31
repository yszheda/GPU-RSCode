// =====================================================================================
// 
//       Filename:  GaloisFieldValue.hpp
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  2013/8/30 ä¸å 11:29:29
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Shuai YUAN (galoisplusplus), yszheda AT gmail DOT com
//        Company:  
// 
// =====================================================================================

#ifndef  GALOISFIELDVALUE_H
#define  GALOISFIELDVALUE_H
#include <iostream>
#include <vector>

enum class GaloisFieldPolicies {LOOP, FULLTABLE, DOUBLETABLES, LOGEXPTABLES, LOGEXPTABLES_V1, LOGEXPTABLES_V2, LOGEXPTABLES_V3};

// =====================================================================================
//        Class:  GaloisFieldValue
//  Description:  
// =====================================================================================

template < unsigned int gf_width = 8 >
class GaloisFieldValue
{
		public:
				// ====================  LIFECYCLE     =======================================
				GaloisFieldValue ();	// default constructor
				GaloisFieldValue ( const int gf_value ); // constructor
				GaloisFieldValue ( const GaloisFieldValue<gf_width> &other ); // copy constructor
				~GaloisFieldValue () { };                          // destructor

				// ====================  ACCESSORS     =======================================

				// ====================  MUTATORS      =======================================

				// ====================  OPERATORS     =======================================

				GaloisFieldValue<gf_width>& operator = ( const GaloisFieldValue<gf_width> &other ); // copy assignment operator
				GaloisFieldValue<gf_width>& operator += ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& operator -= ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& operator *= ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& operator /= ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& operator ^= ( const int &rhs );
				inline bool operator == ( const GaloisFieldValue<gf_width> &other ) const 
				{
						return (this.gf_value == other.gf_value);
				}
				inline bool operator != ( const GaloisFieldValue<gf_width> &other ) const
				{
						return !(*this == other);
				}
				inline bool operator < ( const GaloisFieldValue<gf_width> &other ) const 
				{
						return (this.gf_value < other.gf_value);
				}
				inline bool operator <= ( const GaloisFieldValue<gf_width> &other ) const 
				{
						return ((*this < other) || (*this == other));
				}
				inline bool operator > ( const GaloisFieldValue<gf_width> &other ) const 
				{
						return !(*this <= other);
				}
				inline bool operator >= ( const GaloisFieldValue<gf_width> &other ) const 
				{
						return !(*this < other);
				}

				// NOTE: ostream operator is also a template 
				// gcc complains about the following declaration, absurd?
				// friend std::ostream& operator << <>( std::ostream& os, const GaloisFieldValue<gf_width> &value );
				// use another way like the following declaration:
				// rename the template parameter is necessary,
				// and no default template arguments are allowed
				template < unsigned int width >
				friend std::ostream& operator << ( std::ostream& os, const GaloisFieldValue<width> &value );

		protected:
				// ====================  DATA MEMBERS  =======================================

		private:
				GaloisFieldValue ( GaloisFieldValue<gf_width> &&other );   // move constructor
				GaloisFieldValue<gf_width>& operator = ( GaloisFieldValue<gf_width> &&other ); 		// move assignment operator
				void setupTables ();
				void setupFullTables ();
				void setupDoubleTables ();
				void setupLogExpTables ();
				void setupLogExpTablesV1 ();
				void setupLogExpTablesV2 ();
				void setupLogExpTablesV3 ();

				int loopMul ( const int x, const int y );
				GaloisFieldValue<gf_width>& loopMulAssign ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& fullTableMulAssign ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& doubleTablesMulAssign ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& logExpTablesMulAssign ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& logExpTablesMulAssignV0 ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& logExpTablesMulAssignV1 ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& logExpTablesMulAssignV2 ( const GaloisFieldValue<gf_width> &other );
				GaloisFieldValue<gf_width>& logExpTablesMulAssignV3 ( const GaloisFieldValue<gf_width> &other );

				// ====================  DATA MEMBERS  =======================================
				int gf_value;		// the value of galois field poly
				GaloisFieldPolicies policy;

				// mul/div tables are used for FULLTABLE policy
				static std::vector<int> gf_mul_table;
				static std::vector<int> gf_div_table;
				// used for DOUBLETABLES policy
				static std::vector<int> gf_left_mult_table;
				static std::vector<int> gf_right_mult_table;
				static std::vector<int> gf_left_div_table;
				static std::vector<int> gf_right_div_table;
				// log/exp tables are used for LOGEXPTABLES policy
				static std::vector<int> gf_log_table;
				static std::vector<int> gf_exp_table;

				// constexpr is a C++11 feature
#if __cplusplus > 201100L
				static constexpr unsigned int prim_poly_table[33] = 
#else
				static const unsigned int prim_poly_table[33] = 
#endif
				{ 0, 
				/*  1 */     1, 
				/*  2 */    07,
				/*  3 */    013,
				/*  4 */    023,
				/*  5 */    045,
				/*  6 */    0103,
				/*  7 */    0211,
				/*  8 */    0435,
				/*  9 */    01021,
				/* 10 */    02011,
				/* 11 */    04005,
				/* 12 */    010123,
				/* 13 */    020033,
				/* 14 */    042103,
				/* 15 */    0100003,
				/* 16 */    0210013,
				/* 17 */    0400011,
				/* 18 */    01000201,
				/* 19 */    02000047,
				/* 20 */    04000011,
				/* 21 */    010000005,
				/* 22 */    020000003,
				/* 23 */    040000041,
				/* 24 */    0100000207,
				/* 25 */    0200000011,
				/* 26 */    0400000107,
				/* 27 */    01000000047,
				/* 28 */    02000000011,
				/* 29 */    04000000005,
				/* 30 */    010040000007,
				/* 31 */    020000000011, 
				/* 32 */    00020000007 };  /* Really 40020000007, but we're omitting the high order bit */
								
}; // -----  end of template class GaloisFieldValue  -----

template < unsigned int gf_width = 8 >
GaloisFieldValue<gf_width> operator + ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs );
template < unsigned int gf_width = 8 >
GaloisFieldValue<gf_width> operator -= ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs );
template < unsigned int gf_width = 8 >
GaloisFieldValue<gf_width> operator *= ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs );
template < unsigned int gf_width = 8 >
GaloisFieldValue<gf_width> operator /= ( const GaloisFieldValue<gf_width> &lhs, const GaloisFieldValue<gf_width> &rhs );
template < unsigned int gf_width = 8 >
GaloisFieldValue<gf_width> operator ^= ( const GaloisFieldValue<gf_width> &lhs, const int &rhs );

#endif   // ----- #ifndef GALOISFIELDVALUE_H  -----
