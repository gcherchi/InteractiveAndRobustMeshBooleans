/****************************************************************************
* NFG - Numbers for Geometry                     					        *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2019: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU Lesser General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or (at  *
* your option) any later version.                                           *
*                                                                           *
* This program is distributed in the hope that it will be useful, but       *
* WITHOUT ANY WARRANTY; without even the implied warranty of                *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  *
* General Public License for more details.                                  *
*                                                                           *
* You should have received a copy of the GNU Lesser General Public License  *
* along with this program.  If not, see http://www.gnu.org/licenses/.       *
*                                                                           *
****************************************************************************/

/////////////////////////////////////////////////////////////////////////////
//
// Compilation directives for common compilers/processors
// (minimum compiler version indicated)
// 
// x86-64 gcc 9.3:		-O3 -std=c++2a -frounding-math [-mavx2 -mfma]
// x86-64 clang 12.0.0:	-O3 -std=c++2a -frounding-math [-mavx2 -mfma]
// x86-64 MSVC 19.30:	/O2 /std:c++20 /fp:strict [/arch:AVX2]
// x86-64 ICX 2021 2.0:	-O3 -std=c++20 -fp-model=strict [-mavx2]
// 
// ARM64 gcc 12.1:		-O3 -std=c++2a -frounding-math
// armv8-a clang 16.00:	-O3 -std=c++2a -frounding-math
// arm64 MSVC 19.30:	/O2 /std:c++20 /fp:strict
//
/////////////////////////////////////////////////////////////////////////////

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <fenv.h>
#include <iostream>
#include <climits>
#include <assert.h>
#include <vector>
#include <bitset>
#include <cstdint>
#include <string>
#include <cstring>
#include <algorithm>

#pragma intrinsic(fabs)

//#define USE_GNU_GMP_CLASSES
#ifdef USE_GNU_GMP_CLASSES
#include <gmpxx.h>
#endif

#if INTPTR_MAX == INT64_MAX
#	ifdef __SSE2__
#		define USE_SIMD_INSTRUCTIONS
#	endif
#	ifdef __AVX2__
#		define USE_SIMD_INSTRUCTIONS
#		define USE_AVX2_INSTRUCTIONS
#	endif
#	ifdef SIMDE_ENABLE_NATIVE_ALIASES
#		include <x86/avx2.h>
#		include <x86/fma.h>
#		define USE_SIMD_INSTRUCTIONS
#		define USE_AVX2_INSTRUCTIONS
#	else
#		ifdef USE_SIMD_INSTRUCTIONS
#			ifdef USE_AVX2_INSTRUCTIONS
#				include <immintrin.h>
#			else
#				include <emmintrin.h>
#			endif
#		endif
#	endif
#endif

#ifdef _MSC_VER
#	pragma fenv_access (on)
#   if _MSVC_LANG >= 202002L
#       define STDCPLUSPLUS20
#   endif
#else
#	pragma STDC FENV_ACCESS ON
#	ifdef __AVX2__
#		pragma GCC target("fma")
#	endif
#   if __cplusplus >= 202002L
#       define STDCPLUSPLUS20
#   endif
#endif

inline void setFPUModeToRoundUP() { fesetround(FE_UPWARD); }
inline void setFPUModeToRoundNEAR() { fesetround(FE_TONEAREST); }

// Deprecated. This function exists for back-compatibility only.
inline void initFPU() { }

inline void ip_error(const char* msg)
{
	std::cerr << msg;
	exit(0);
}

#ifdef STDCPLUSPLUS20
#include <bit>
inline int nfg_count_lz(uint32_t v) { return std::countl_zero(v); }
inline int nfg_count_rz(uint32_t v) { return std::countr_zero(v); }
#else
// Slower versions of the above functions for C++ standards < 20
inline int nfg_count_lz(uint32_t v) {
	int z = 32;
	while (v) { v >>= 1; z--; }
	return z;
}

inline int nfg_count_rz(uint32_t v) {
	int z = 32;
	while (v) { v <<= 1; z--; }
	return z;
}
#endif

/////////////////////////////////////////////////////////////////////
// 	   
// 	   M E M O R Y   P O O L S
// 
/////////////////////////////////////////////////////////////////////

// An N_block is a contiguous portion of memory elements.
// An N_memory_pool stores a set of N_blocks along with
// a 'stack' of pointers to these blocks.
// Upon allocation, the first free block is picked from the stack.
// If the stack is empty, the size of the pool (and the stack) is grown.
// When a block is released, its pointer is added on top of the stack.
//
// alloc is almost always O(1) [O(num_blocks) only if growth is required]
// release is always O(1)

class N_memory_pool {
	std::vector<uint8_t*> data; // Data storage. Each data[i] is an array of N_blocks
	uint8_t** stack;	// Stack of pointers to free blocks
	size_t last;		// Element past to the last in the stack
	size_t size;		// Size of the pool (total number of blocks)
	const size_t block_size; // Size of one block (number N of elements)

public:
	N_memory_pool(size_t _size, size_t _block_size);
	N_memory_pool(N_memory_pool&& p) noexcept;
	~N_memory_pool();

	void* alloc();
	void release(void* s);

protected:
	void addDataBlockArray(size_t tot_size, uint8_t** fs);
	void doubleDataArray();
};

// An N_block is a contiguous portion of 32-bit memory elements.
// An extended_N_block is the concatentaion of [N] and an N_block.
//
//		i.e. uint32_t block_N[N+1] = { N, v1, v2, ..., vN };
// 
// An Indexed N pool (IN_pool) stores a set of extended_N_blocks along with
// a 'stack' of pointers to these extended blocks. Each pointer in the stack
// points to the second element of an extended block because the first must 
// always store the block size and is therefore reserved.
// Upon allocation, the first free block is picked from the stack.
// If the stack is empty, the size of the pool (and the stack) is grown.
// When a block is released, its pointer is added on top of the stack.
//
// alloc is almost always O(1) [O(num_blocks) only if growth is required]
// release is always O(1)

class IN_pool {
	std::vector<uint32_t*> data; // Data storage. Each data[i] is an array of extended_N_blocks
	uint32_t** stack;	// Stack of pointers to free blocks
	uint32_t last;		// Element past to the last in the stack
	uint32_t size;		// Size of the pool (total number of blocks)
	const uint32_t block_size; // Size of one block (number N of elements)

public:
	IN_pool(uint32_t _size, uint32_t _block_size);
	IN_pool(IN_pool&& p) noexcept;
	~IN_pool();

	uint32_t* alloc();
	void release(uint32_t* s);

	uint32_t blockSize() const;

protected:
	void addDataBlockArray(uint32_t tot_size, uint32_t extended_block_size, uint32_t** fs);
	void doubleDataArray();
};

// A MultiPool is a collection of IN_pools having size N = 2, 4, 8, 16, ..., 2^m
// 
// Upon allocation for X elements two cases may occur:
// 1) X <= 2^m
//     Calculate the minimum value of N such that N >= X
//     Allocate a block in the corresponding IN_pool
// 2) X > 2^m 
//	   standard malloc() is used for X+1 elements
// 	   first element is set to 0 [meaning that no IN_pool was used]
// 	   a pointer to the second element is returned for use
//
// Upon release of a pointer A
//    Let V be the value of the element preceeding the one pointed by A (i.e. V = A[-1])
//    If V==0 use standard free()
//    else if V==N release a block is the corresponding IN_pool

class MultiPool {
	std::vector<IN_pool> IN_pools;
	uint32_t max_block_size;

public:
	MultiPool(uint32_t _max_block_size = 32, uint32_t init_capacity = 16);

	void* alloc(uint32_t num_bytes);
	void release(void* s);

protected:
	IN_pool& pickPoolFromSize(uint32_t bs);
};


/////////////////////////////////////////////////////////////////////
// 	   
// 	   I N T E R V A L   A R I T H M E T I C
// 
/////////////////////////////////////////////////////////////////////

// An interval_number is a pair of doubles representing an interval.
// Operations on interval_number require that the rounding mode is
// set to +INFINITY. Use setFPUModeToRoundUP().

class interval_number
{
#ifdef USE_SIMD_INSTRUCTIONS
	__m128d interval; // interval[1] = min_low, interval[0] = high
#else // USE_SIMD_INSTRUCTIONS
	double min_low, high;
#endif // USE_SIMD_INSTRUCTIONS

public:
	interval_number();
	interval_number(const double a);
	interval_number(const double minf, const double sup);
	interval_number(const interval_number& b);

	const double* getInterval() const;

	double minus_inf() const;
	double inf() const;
	double sup() const;

	bool isNegative() const;
	bool isPositive() const;
	void negate();

	interval_number& operator=(const interval_number& b);

	int sign() const; // Zero is not accounted for

	interval_number operator-() const;

	interval_number operator+(const interval_number& b) const;
	interval_number operator-(const interval_number& b) const;
	interval_number operator*(const interval_number& b) const;
	interval_number operator+(const double b) const;
	interval_number operator-(const double b) const;
	interval_number operator*(const double b) const;
	interval_number operator/(const double b) const;
	interval_number& operator+=(const interval_number& b);
	interval_number& operator-=(const interval_number& b);
	interval_number& operator*=(const interval_number& b);
	interval_number& operator+=(const double b);
	interval_number& operator-=(const double b);
	interval_number& operator*=(const double b);
	interval_number& operator/=(const double b);

	interval_number abs() const;
	interval_number sqr() const;
	interval_number pow(unsigned int e) const;
	friend interval_number min(const interval_number& a, const interval_number& b);
	friend interval_number max(const interval_number& a, const interval_number& b);

	double width() const;

	bool signIsReliable() const; // Zero is not accounted for
	bool containsZero() const;

	bool isNAN() const;

	double getMid() const;
	bool isExact() const;

	// Can be TRUE only if the intervals are disjoint
	bool operator<(const interval_number& b) const;
	bool operator>(const interval_number& b) const;

	// Can be TRUE only if the interval interiors are disjoint
	bool operator<=(const interval_number& b) const;
	bool operator>=(const interval_number& b) const;

	// TRUE if the intervals are identical single values
	bool operator==(const interval_number& b) const;

	// TRUE if the intervals have no common values
	bool operator!=(const interval_number& b) const;

	bool operator<(const double b) const;
	bool operator<=(const double b) const;
	bool operator>(const double b) const;
	bool operator>=(const double b) const;
	bool operator==(const double b) const;

	// The inverse of an interval. Returns NAN if the interval contains zero
	interval_number inverse() const;

#ifdef USE_SIMD_INSTRUCTIONS
	interval_number(const __m128d& i);

protected:
	static __m128d zero();
	static __m128d minus_one();
	static __m128d sign_low_mask();
	static __m128d sign_high_mask();
	static __m128d sign_fabs_mask();
	static __m128d all_high_mask();

	__m128d getLowSwitched() const;
#endif // USE_SIMD_INSTRUCTIONS
};

// The square root of an interval
// Returns NAN if the interval contains a negative value
interval_number sqrt(const interval_number& p);

// stream output
inline std::ostream& operator<<(std::ostream& os, const interval_number& p)
{
	os << "[ " << p.inf() << ", " << p.sup() << " ]";
	return os;
}


/////////////////////////////////////////////////////////////////////
// 	   
// 	   E X P A N S I O N   A R I T H M E T I C
// 
/////////////////////////////////////////////////////////////////////

// An instance of the following must be created to access functions for expansion arithmetic
struct expansionObject
{
public:
	inline static thread_local MultiPool mempool = MultiPool(2048, 64);

	// Basic error-free transformations
	static void Quick_Two_Sum(const double a, const double b, double& x, double& y);
	static void quick_Two_Sum(const double a, const double b, double* xy);
	static void Two_Sum(const double a, const double b, double& x, double& y);
	static void two_Sum(const double a, const double b, double* xy);
	static void Two_One_Sum(const double a1, const double a0, const double b, double& x2, double& x1, double& x0);
	static void two_One_Sum(const double *a, const double b, double* x);
	static void Two_Diff(const double a, const double b, double& x, double& y);
	static void two_Diff(const double a, const double b, double* xy);
	static void Two_One_Diff(const double a1, const double a0, const double b, double& x2, double& x1, double& x0);

	// Products

	// [x,y] = [a]*[b]		 Multiplies two expansions [a] and [b] of length one
	static void Two_Prod(const double a, const double b, double& x, double& y);
	static void Two_Prod(const double a, const double b, double* xy);

	// [x,y] = [a]^2		Squares an expansion of length one
	static void Square(const double a, double& x, double& y);
	static void Square(const double a, double* xy);

	// [x2,x1,x0] = [a1,a0]-[b]		Subtracts an expansion [b] of length one from an expansion [a1,a0] of length two
	static void two_One_Diff(const double a1, const double a0, const double b, double& x2, double& x1, double& x0);
	static void two_One_Diff(const double* a, const double b, double* x);

	// [x3,x2,x1,x0] = [a1,a0]*[b]		Multiplies an expansion [a1,a0] of length two by an expansion [b] of length one
	static void Two_One_Prod(const double a1, const double a0, const double b, double& x3, double& x2, double& x1, double& x0);
	static void Two_One_Prod(const double* a, const double b, double* x);

	// [x3,x2,x1,x0] = [a1,a0]+[b1,b0]		Calculates the sum of two expansions of length two
	static void Two_Two_Sum(const double a1, const double a0, const double b1, const double b0, double& x3, double& x2, double& x1, double& x0);
	static void Two_Two_Sum(const double* a, const double* b, double* xy);

	// [x3,x2,x1,x0] = [a1,a0]-[b1,b0]		Calculates the difference between two expansions of length two
	static void Two_Two_Diff(const double a1, const double a0, const double b1, const double b0, double& x3, double& x2, double& x1, double& x0);
	static void Two_Two_Diff(const double* a, const double* b, double* x);

	// Calculates the second component 'y' of the expansion [x,y] = [a]-[b] when 'x' is known
	static void Two_Diff_Back(const double a, const double b, double& x, double& y);
	static void Two_Diff_Back(const double a, const double b, double* xy);

	// [h] = [a1,a0]^2		Squares an expansion of length 2
	// 'h' must be allocated by the caller with 6 components.
	static void Two_Square(const double& a1, const double& a0, double* x);
	static void Two_Square(const double* a, double* x);

	// [h7,h6,...,h0] = [a1,a0]*[b1,b0]		Calculates the product of two expansions of length two.
	// 'h' must be allocated by the caller with eight components.
	static void Two_Two_Prod(const double a1, const double a0, const double b1, const double b0, double* h);
	static void Two_Two_Prod(const double* a, const double* b, double* xy);

	// [e] = -[e]		Inplace inversion
	static void Gen_Invert(const int elen, double* e);

	// [h] = [e] + [f]		Sums two expansions and returns number of components of result
	// 'h' must be allocated by the caller with at least elen+flen components.
	static int Gen_Sum(const int elen, const double* e, const int flen, const double* f, double* h);

	// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
	static int Gen_Sum_With_Alloc(const int elen, const double* e, const int flen, const double* f, double** h);

	// [h] = [e] + [f]		Subtracts two expansions and returns number of components of result
	// 'h' must be allocated by the caller with at least elen+flen components.
	static int Gen_Diff(const int elen, const double* e, const int flen, const double* f, double* h);

	// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
	static int Gen_Diff_With_Alloc(const int elen, const double* e, const int flen, const double* f, double** h);

	// [h] = [e] * b		Multiplies an expansion by a scalar
	// 'h' must be allocated by the caller with at least elen*2 components.
	static int Gen_Scale(const int elen, const double* e, const double b, double* h);

	// [h] = [e] * 2		Multiplies an expansion by 2
	// 'h' must be allocated by the caller with at least elen components. This is exact up to overflows.
	static void Double(const int elen, const double* e, double* h);

	// [h] = [e] * n		Multiplies an expansion by n
	// If 'n' is a power of two, the multiplication is exact
	static void ExactScale(const int elen, double* e, const double n);

	// [h] = [a] * [b]
	// 'h' must be allocated by the caller with at least 2*alen*blen components.
	static int Sub_product(const int alen, const double* a, const int blen, const double* b, double* h);

	// [h] = [a] * [b]
	// 'h' must be allocated by the caller with at least MAX(2*alen*blen, 8) components.
	static int Gen_Product(const int alen, const double* a, const int blen, const double* b, double* h);

	// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
	static int Gen_Product_With_Alloc(const int alen, const double* a, const int blen, const double* b, double** h);


	// Assume that *h is pre-allocated with hlen doubles.
	// If more elements are required, *h is re-allocated internally.
	// In any case, the function returns the size of the resulting expansion.
	// The caller must verify whether reallocation took place, and possibly call 'free' to release the memory.
	// When reallocation takes place, *h becomes different from its original value.
	static int Double_With_PreAlloc(const int elen, const double* e, double** h, const int hlen);
	static int Gen_Scale_With_PreAlloc(const int elen, const double* e, const double& b, double** h, const int hlen);
	static int Gen_Sum_With_PreAlloc(const int elen, const double* e, const int flen, const double* f, double** h, const int hlen);
	static int Gen_Diff_With_PreAlloc(const int elen, const double* e, const int flen, const double* f, double** h, const int hlen);
	static int Gen_Product_With_PreAlloc(const int alen, const double* a, const int blen, const double* b, double** h, const int hlen);

	// Approximates the expansion to a double
	static double To_Double(const int elen, const double* e);

	static void print(const int elen, const double* e);

#ifndef USE_AVX2_INSTRUCTIONS 
	static void Split(double a, double& _ah, double& _al);
	static void Two_Prod_PreSplit(double a, double b, double _bh, double _bl, double& x, double& y);
	static void Two_Product_2Presplit(double a, double _ah, double _al, double b, double _bh, double _bl, double& x, double& y);
#endif
};

#ifdef USE_GNU_GMP_CLASSES
typedef mpz_class bignatural;
typedef mpq_class bigfloat;
typedef mpq_class bigrational;

inline void read_rational(FILE* fp, bigrational& r) {
	gmp_fscanf(fp, "%Qd,", r.get_mpq_t()); ungetc(',', fp);
}

inline bigfloat getBigFloatFromRational(const bigrational& r, uint32_t prec_bits) { return r; }

inline bigfloat sqrt(const bigfloat& f, uint32_t prec_bits) {
	mpf_class bf(f, prec_bits);
	mpf_class s = sqrt(bf);
	return mpq_class(s);
}

inline int32_t log2(const bigfloat& f) {
	mpf_class bf(f);
	return bf.get_mpf_t()->_mp_exp + (int32_t)bf.get_prec() - 1;
}

inline void add1ULP(bigfloat& f) {
	mpf_class bf(f);
	mpf_class ulp(1U, bf.get_prec());
	ulp.get_mpf_t()->_mp_exp = bf.get_mpf_t()->_mp_exp;
	bf += ulp;
}

#else
/////////////////////////////////////////////////////////////////////
// 	   
// 	   B I G   N A T U R A L
// 
/////////////////////////////////////////////////////////////////////

// A bignatural is an arbitrarily large non-negative integer.
// It is made of a sequence of digits in base 2^32.
// Leading zero-digits are not allowed.
// The value 'zero' is represented by an empty digit sequence.

class bignatural {
	uint32_t m_capacity;	// Current vector capacity
	uint32_t m_size;		// Actual number of digits
	uint32_t* digits;	    // Ptr to the digits

	static uint32_t* BN_ALLOC(uint32_t num_bytes);
	static void BN_FREE(uint32_t* ptr);

	// Read as many decimal digits as possible from file so that they fit a uint64_t
	// Return the number of digits read
	static size_t scan_uint64_t(FILE* fp, uint64_t& t);

	void init(const bignatural& m);
	void init(const uint32_t m);
	void init(const uint64_t m);
	void init(FILE* fp);

public:
	// Creates a 'zero'
	bignatural();

	// Destructor
	~bignatural();

	// Copy constructor
	bignatural(const bignatural& m);

	// Move constructor
	bignatural(bignatural&& m) noexcept;

	// Construct from unsigned 32bit integer
	bignatural(uint32_t m);

	// Construct from unsigned 64bit integer
	bignatural(uint64_t m);

	// Construct from FILE stream
	bignatural(FILE* f);

	// If the number fits a uint64_t convert and return true
	bool toUint64(uint64_t& n) const;

	// If the number fits a uint32_t convert and return true
	bool toUint32(uint32_t& n) const;

	// Assignment operators
	bignatural& operator=(const bignatural& m);
	bignatural& operator=(const uint64_t m);

	// Get the least significant digit
	const uint32_t& back() const;

	// Get the i'th digit
	const uint32_t& operator[](int i) const;

	// Number of significant digits
	uint32_t size() const;

	// TRUE if number is zero
	bool empty() const;

	// Number of significant bits
	uint32_t getNumSignificantBits() const;

	// Get the i'th bit
	bool getBit(uint32_t i) const;

	// Left-shift by n bits and possibly add limbs as necessary
	void operator<<=(uint32_t n);
	bignatural operator<<(uint32_t n) const;

	// Right-shift by n bits
	void operator>>=(uint32_t n);
	bignatural operator>>(uint32_t n) const;

	// Comparison operators
	bool operator==(const bignatural& b) const;
	bool operator!=(const bignatural& b) const;
	bool operator>=(const bignatural& b) const;
	bool operator>(const bignatural& b) const;
	bool operator<=(const bignatural& b) const;
	bool operator<(const bignatural& b) const;

	// Arithmetic operations
	bignatural& operator+=(const bignatural& b);
	bignatural& operator+=(const uint32_t b);
	bignatural& operator+=(const uint64_t b);
	bignatural operator+(const bignatural& b) const;
	bignatural operator+(const uint32_t b) const;
	bignatural operator+(const uint64_t b) const;

	// Assume that b is smaller than or equal to this number!
	bignatural& operator-=(const bignatural& b);
	bignatural operator-(const bignatural& b) const;

	bignatural operator*(const bignatural& b) const;
	bignatural& operator*=(const bignatural& b);
	bignatural& operator*=(const uint32_t b);
	bignatural& operator*=(const uint64_t b);

	// Short division
	bignatural divide_by(const uint32_t D, uint32_t& remainder) const;

	// Long division
	bignatural divide_by(const bignatural& divisor, bignatural& remainder) const;

	// Long sqrt (truncated)
	bignatural sqrt() const;

	// Bitwise OR
	bignatural operator|(const bignatural& b) const;
	void operator|=(uint32_t i);

	// Greatest common divisor (Euclidean algorithm)
	bignatural GCD(const bignatural& D) const;

	// String representation in decimal form
	std::string get_dec_str() const;

	// String representation in binary form
	std::string get_str() const;

	// Count number of zeroes on the right (least significant binary digits)
	uint32_t countEndingZeroes() const;

protected:
	uint32_t& back();

	void pop_back();

	uint32_t& operator[](int i);

	void push_back(uint32_t b);

	// Left-shift. Same as above but assumes that number is not zero!
	void leftShift(uint32_t n);

	void push_bit_back(uint32_t b);

	void reserve(uint32_t n);

	void resize(uint32_t n);

	void fill(uint32_t v);

	void pop_front();

	// Count number of zeroes on the left (most significant digits).
	// Assumes that number is not zero!
	uint32_t countLeadingZeroes() const;

	// Count number of zeroes on the right of the last 1 in the least significant limb
	// Assumes that number is not zero and last limb is not zero!
	uint32_t countEndingZeroesLSL() const;

	void pack();

	// a and b must NOT be this number!
	void toSum(const bignatural& a, const bignatural& b);

	// a and b must NOT be this number!
	// Assume that b is smaller or equal than a!
	void toDiff(const bignatural& a, const bignatural& b);

	// a and b must NOT be this number!
	void toProd(const bignatural& a, const bignatural& b);

private:

	// Multiplies by a single limb, left shift, and add to accumulator. Does not pack!
	void addmul(uint32_t b, uint32_t left_shifts, bignatural& result) const;

	// Increases the vector capacity while maintaining the number validity
	void increaseCapacity(uint32_t new_capacity);

	// Adds one most significant digit while making room if necessary
	void addOneMostSignificantDigit(uint32_t d);

	// Memory pool for bignaturals.
	inline static thread_local MultiPool nfgMemoryPool;

	friend class bigfloat;
	friend class bigrational;
};

// Operators with left-hand primitive integers
inline bignatural operator+(uint32_t a, const bignatural& p) { return p + a; }
inline bignatural operator+(uint64_t a, const bignatural& p) { return p + a; }
inline bignatural operator*(uint64_t a, const bignatural& p) { return p * a; }

inline bignatural sqrt(const bignatural& n) { return n.sqrt(); }

inline std::ostream& operator<<(std::ostream& os, const bignatural& p)
{
	os << p.get_dec_str();
	return os;
}

/////////////////////////////////////////////////////////////////////
// 	   
// 	   B I G   F L O A T
// 
/////////////////////////////////////////////////////////////////////

// A bigfloat is a floting point number with arbitrarily large mantissa.
// In principle, we could have made the exponent arbitrarily large too,
// but in practice this appears to be useless.
// Exponents are in the range [-INT32_MAX, INT32_MAX]
//
// A bigfloat f evaluates to f = sign * mantissa * 2^exponent
//
// mantissa is a bignatural whose least significant bit is 1.
// Number is zero if mantissa is empty.

class bigfloat {
	bignatural mantissa; // .back() is less significant. Use 32-bit limbs to avoid overflows using 64-bits
	int32_t exponent; // In principle we might still have under/overflows, but not in practice
	int32_t sign;	// Redundant but keeps alignment

public:
	// Default constructor creates a zero-valued bigfloat
	bigfloat();

	// Lossless conversion from double
	bigfloat(const double d);

	// Constructs from a bignatural
	bigfloat(const bignatural& m, int32_t e, int32_t s);

	// Truncated approximation
	double get_d() const;

	// Compute a sqrt as precise as prec_bits bits (num. bits in resulting mantissa)
	bigfloat sqrt(uint32_t prec_bits) const;

	// Truncated base-2 logarithm
	int32_t log2() const;

	// Adds one ULP to the number
	void increaseMantissa();

	// Arithmetic operations
	bigfloat operator+(const bigfloat& b) const;
	bigfloat operator-(const bigfloat& b) const;
	bigfloat operator*(const bigfloat& b) const;

	// Comparison
	bool operator!=(const bigfloat& b) const;

	// Sign switch
	void invert();

	bigfloat inverse() const;

	// Get sign
	int sgn() const;

	// Convert to string (binary exponential representation)
	std::string get_str() const;

	// Access components
	const bignatural& getMantissa() const;
	int32_t getExponent() const;

private:
	// Right-shift as long as the least significant bit is zero
	void pack();

	// Left-shift the mantissa by n bits and reduce the exponent accordingly
	void leftShift(uint32_t n);

	// Right-shift the mantissa by n bits and reduce the exponent accordingly
	void rightShift(uint32_t n);
};

// Sign operator for bigfloats
inline int sgn(const bigfloat& f) { return f.sgn(); }

// Operators with left-hand doubles
inline bigfloat operator+(const double a, const bigfloat& p) { return p + a; }
inline bigfloat operator-(const double a, const bigfloat& p) { return (p - a).inverse(); }
inline bigfloat operator*(const double a, const bigfloat& p) { return p * a; }

// Sign inversion operator for bigfloats
inline bigfloat operator-(const bigfloat& f) { return f.inverse(); }

inline bigfloat sqrt(const bigfloat& f, uint32_t prec_bits) {
	return f.sqrt(prec_bits);
}

inline int32_t log2(const bigfloat& f) { return f.log2(); }

inline void add1ULP(bigfloat& f) { f.increaseMantissa(); }

// std::ostream interface (decimal exponential representation)
inline std::ostream& operator<<(std::ostream& os, const bigfloat& p) {
	if (p.sgn() < 0) os << "-";
	os << p.getMantissa().get_dec_str() << " * 2^" << p.getExponent();
	return os;
}


/////////////////////////////////////////////////////////////////////
// 	   
// 	   B I G   R A T I O N A L
// 
/////////////////////////////////////////////////////////////////////

// A bigrational is a fraction of two bignaturals with a sign.
// Number is zero if sign is zero

class bigrational {
	bignatural numerator, denominator;
	int32_t sign;

public:
	// Create a zero
	bigrational() : sign(0) {}

	// Create from a double (lossless)
	bigrational(const double f) : bigrational(bigfloat(f)) {}

	// Create from a bigfloat (lossless)
	bigrational(const bigfloat& f);

	// Create from a file
	bigrational(FILE* fp) { init(fp); }

	// Create from explicit numerator, denominator and sign.
	bigrational(const bignatural& num, const bignatural& den, int32_t s) :
		numerator(num), denominator(den), sign(s) {
		canonicalize();
	}

	// Convert to multiplicative inverse
	void invert();

	// Return multiplicative inverse
	bigrational inverse() const;

	// Invert sign
	void negate();

	// Return additive inverse
	bigrational negation() const;

	// Arithmetic operations
	bigrational operator+(const bigrational& r) const;
	bigrational operator-(const bigrational& r) const;
	bigrational operator*(const bigrational& r) const;
	bigrational operator/(const bigrational& r) const;

	// Comparison operators
	bool operator==(const bigrational& r) const;
	bool operator!=(const bigrational& r) const;
	bool operator>(const bigrational& r) const;
	bool operator>=(const bigrational& r) const;
	bool operator<(const bigrational& r) const;
	bool operator<=(const bigrational& r) const;
	bool hasGreaterModule(const bigrational& r) const;
	bool hasGrtrOrEqModule(const bigrational& r) const;

	// Conversion to double (truncated)
	double get_d() const;

	// Conversion to bigfloat (truncated)
	bigfloat get_bigfloat(uint32_t num_significant_bits) const;

	// Access to components
	const bignatural& get_num() const;
	const bignatural& get_den() const;

	// Get sign
	int32_t sgn() const;

	// Return decimal representation
	std::string get_dec_str() const;

	// Return binary representation
	std::string get_str() const;

protected:
	// Iteratively divide both num and den by two as long as they are both even
	void compress();

	// Make numerator and denominator coprime (divide both by GCD)
	void canonicalize();

	void init(FILE* fp);
};

// Operators with left-hand doubles
inline bigrational operator+(const double a, const bigrational& p) { return p + a; }
inline bigrational operator-(const double a, const bigrational& p) { return (p - a).negation(); }
inline bigrational operator*(const double a, const bigrational& p) { return p * a; }
inline bigrational operator/(const double a, const bigrational& p) { return bigrational(a) / p; }

// Sign inversion operator for bigrationals
inline bigrational operator-(const bigrational& p) { return p.negation(); }

// Sign operator for bigrationals
inline int32_t sgn(const bigrational& p) { return p.sgn(); }

// std::ostream interface (decimal fraction representation)
inline std::ostream& operator<<(std::ostream& os, const bigrational& p)
{
	os << p.get_dec_str();
	return os;
}

inline void read_rational(FILE* fp, bigrational& r) {
	r = bigrational(fp);
}

inline bigfloat getBigFloatFromRational(const bigrational& r, uint32_t prec_bits) { return r.get_bigfloat(prec_bits); }

#endif // USE_GNU_GMP_CLASSES


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// FUNCTION IMPLEMENTATIONS - MEMORY POOLS
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

inline void N_memory_pool::addDataBlockArray(size_t tot_size, uint8_t** fs) {
	uint8_t* data_1 = (uint8_t*)malloc(sizeof(uint8_t) * tot_size);
	data.push_back(data_1);
	uint8_t* lst = data_1 + tot_size;
	while (lst != data_1) {
		lst -= block_size;
		*fs++ = lst;
	}
}

inline void N_memory_pool::doubleDataArray() {
	uint8_t** fs = (uint8_t**)malloc(sizeof(uint8_t*) * size * 2);
	size_t i = 0, j = size;
	while (i < size) fs[j++] = stack[i++];
	free(stack);
	stack = fs;

	addDataBlockArray(block_size * size, fs);

	last += size;
	size *= 2;
}

inline N_memory_pool::N_memory_pool(size_t _size, size_t _block_size) : last(_size), size(_size), block_size(_block_size) {
	stack = (uint8_t**)malloc(sizeof(uint8_t*) * size);
	uint8_t** fs = stack;
	addDataBlockArray(size * block_size, fs);
}

inline N_memory_pool::~N_memory_pool() {
	free(stack);
	for (uint8_t* p : data) free(p);
	data.clear();
}

inline void* N_memory_pool::alloc() {
	if (last == 0) doubleDataArray();
	return stack[--last];
}

inline void N_memory_pool::release(void* s) {
	stack[last++] = (uint8_t*)s;
}

inline N_memory_pool::N_memory_pool(N_memory_pool&& p) noexcept
	: data(p.data), stack(p.stack), last(p.last), size(p.size), block_size(p.block_size) {
	p.stack = nullptr;
	p.data.clear();
}


inline void IN_pool::addDataBlockArray(uint32_t tot_size, uint32_t extended_block_size, uint32_t** fs) {
	uint32_t* data_1 = (uint32_t*)malloc(sizeof(uint32_t) * tot_size);
	data.push_back(data_1);
	uint32_t* lst = data_1 + tot_size;
	while (lst != data_1) {
		lst -= extended_block_size;
		*lst = block_size;
		*fs++ = lst + 1;
	}
}

inline void IN_pool::doubleDataArray() {
	const uint32_t extended_block_size = block_size + 1;

	uint32_t** fs = (uint32_t**)malloc(sizeof(uint32_t*) * size * 2);
	for (uint32_t i = 0; i < size; i++) fs[i + size] = stack[i];
	free(stack);
	stack = fs;

	addDataBlockArray(extended_block_size * size, extended_block_size, fs);

	last += size;
	size *= 2;
}

inline IN_pool::IN_pool(uint32_t _size, uint32_t _block_size) : last(_size), size(_size), block_size(_block_size) {
	const uint32_t extended_block_size = block_size + 1;
	stack = (uint32_t**)malloc(sizeof(uint32_t*) * size);
	uint32_t** fs = stack;
	addDataBlockArray(size * extended_block_size, extended_block_size, fs);
}


inline IN_pool::~IN_pool() {
	free(stack);
	for (uint32_t* p : data) free(p);
	data.clear();
}

inline uint32_t* IN_pool::alloc() {
	if (last == 0) doubleDataArray();
	return stack[--last];
}

inline void IN_pool::release(uint32_t* s) {
	stack[last++] = s;
}

inline uint32_t IN_pool::blockSize() const {
	return block_size;
}

inline IN_pool::IN_pool(IN_pool&& p) noexcept
	: data(p.data), stack(p.stack), last(p.last), size(p.size), block_size(p.block_size) {
	p.stack = nullptr;
	p.data.clear();
}


inline IN_pool& MultiPool::pickPoolFromSize(uint32_t bs) {
	uint32_t cz = (uint32_t)nfg_count_lz(bs - 1);
	cz = (cz == 32) ? (31) : (cz);
	return IN_pools[31U - cz];
}

inline MultiPool::MultiPool(uint32_t _max_block_size, uint32_t init_capacity) : max_block_size(_max_block_size) {
	uint32_t bs = 1;
	while (bs < max_block_size) {
		bs <<= 1;
		IN_pools.push_back(IN_pool(init_capacity, bs));
	}
}

inline void* MultiPool::alloc(uint32_t num_bytes) {
	const uint32_t num_els = ((num_bytes + 3) >> 2);
	if (num_els > max_block_size) {
		uint32_t* ptr;
		if ((ptr = (uint32_t*)malloc(sizeof(uint32_t) * (num_els + 1))) == NULL) return NULL;
		ptr[0] = 0;
		return ptr + 1;
	}

	return pickPoolFromSize(num_els).alloc();
}

inline void MultiPool::release(void* s) {
	if (s) {
		uint32_t* ps = (((uint32_t*)s) - 1);
		const uint32_t v = *ps;
		if (v == 0) {
			free(ps);
			return;
		}

		pickPoolFromSize(v).release(((uint32_t*)s));
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// FUNCTION IMPLEMENTATIONS - INTERVAL NUMBER
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_SIMD_INSTRUCTIONS
inline __m128d interval_number::zero() { return _mm_setzero_pd(); }
inline __m128d interval_number::minus_one() { return _mm_set1_pd(-1.0); }
inline __m128d interval_number::sign_low_mask() { return _mm_castsi128_pd(_mm_set_epi64x(LLONG_MIN, 0)); }
inline __m128d interval_number::sign_high_mask() { return _mm_castsi128_pd(_mm_set_epi64x(0, LLONG_MIN)); }
inline __m128d interval_number::sign_fabs_mask() { return _mm_castsi128_pd(_mm_set_epi64x(~LLONG_MIN, ~LLONG_MIN)); }
inline __m128d interval_number::all_high_mask() { return _mm_castsi128_pd(_mm_set_epi64x(0, -1LL)); }

inline __m128d interval_number::getLowSwitched() const { return _mm_xor_pd(interval, sign_low_mask()); }

inline interval_number::interval_number() { }
inline interval_number::interval_number(const double a) : interval(_mm_set_pd(-a, a)) { }
inline interval_number::interval_number(const double minf, const double sup) : interval(_mm_set_pd(minf, sup)) { }
inline interval_number::interval_number(const __m128d& i) : interval(i) { }
inline interval_number::interval_number(const interval_number& b) : interval(b.interval) { }

inline double interval_number::minus_inf() const { return _mm_cvtsd_f64(_mm_shuffle_pd(interval, interval, 1)); }
inline double interval_number::inf() const { return -minus_inf(); }
inline double interval_number::sup() const { return _mm_cvtsd_f64(interval); }

inline const double* interval_number::getInterval() const { return (const double*)&interval; }
inline interval_number& interval_number::operator=(const interval_number& b) { interval = b.interval; return *this; }
inline interval_number interval_number::operator+(const interval_number& b) const { return interval_number(_mm_add_pd(interval, b.interval)); }
inline interval_number interval_number::operator-(const interval_number& b) const { return interval_number(_mm_add_pd(interval, _mm_shuffle_pd(b.interval, b.interval, 1))); }
inline interval_number interval_number::operator-() const { return interval_number(_mm_shuffle_pd(interval, interval, 1)); }
inline interval_number interval_number::operator+(const double b) const { return interval_number(_mm_add_pd(interval, _mm_set_pd(-b, b))); }
inline interval_number interval_number::operator-(const double b) const { return interval_number(_mm_sub_pd(interval, _mm_set_pd(-b, b))); }
inline interval_number& interval_number::operator+=(const interval_number& b) { return operator=(*this + b); }
inline interval_number& interval_number::operator-=(const interval_number& b) { return operator=(*this - b); }
inline interval_number& interval_number::operator*=(const interval_number& b) { return operator=(*this * b); }
inline interval_number& interval_number::operator+=(const double b) { return operator=(*this + b); }
inline interval_number& interval_number::operator-=(const double b) { return operator=(*this - b); }
inline interval_number& interval_number::operator*=(const double b) { return operator=(*this * b); }
inline interval_number& interval_number::operator/=(const double b) { return operator=(*this / b); }
inline bool interval_number::operator<(const double b) const { return sup() < b; }
inline bool interval_number::operator<=(const double b) const { return sup() <= b; }
inline bool interval_number::operator>(const double b) const { return inf() > b; }
inline bool interval_number::operator>=(const double b) const { return inf() >= b; }
inline bool interval_number::operator==(const double b) const { return sup() == inf() && sup() == b; }

inline void interval_number::negate() { interval = _mm_shuffle_pd(interval, interval, 1); }

inline bool interval_number::isNegative() const { return _mm_comilt_sd(interval, zero()); }
inline bool interval_number::isPositive() const { return _mm_comilt_sd(_mm_shuffle_pd(interval, interval, 1), zero()); }
inline int interval_number::sign() const { return (isNegative()) ? (-1) : (1); } // Zero is not accounted for


#ifdef USE_AVX2_INSTRUCTIONS
inline interval_number interval_number::operator*(const interval_number& b) const
{
	// This version exploits 256bit registers provided by AVX2 architectures
	// to compute the product using the "naive" eight-multiplications method.
	// The advantage wrt to the non-avx version is due to the absence of
	// branches in the execution, which increses the processor's throughput.

	// Fill i1 and i2 with two copies of 'this' and 'b' respectively
	__m256d i1 = _mm256_castpd128_pd256(interval);
	i1 = _mm256_insertf128_pd(i1, interval, 1);
	__m256d i2 = _mm256_castpd128_pd256(b.interval);
	i2 = _mm256_insertf128_pd(i2, b.interval, 1);

	// Swizzle and change sign appropriately to produce all the eight configs
	__m256d x2 = _mm256_shuffle_pd(i1, i1, 5);
	__m256d x3 = _mm256_xor_pd(i2, _mm256_set_pd(-0.0, -0.0, 0.0, 0.0));
	__m256d x4 = _mm256_xor_pd(i2, _mm256_set_pd(0.0, 0.0, -0.0, -0.0));
	x3 = _mm256_mul_pd(i1, x3);
	x2 = _mm256_mul_pd(x2, x4);
	x3 = _mm256_max_pd(x3, x2);
	x4 = _mm256_shuffle_pd(x3, x3, 5);
	x3 = _mm256_max_pd(x3, x4);
	x3 = _mm256_permute4x64_pd(x3, 72);

	// The first two vals of the 256 vector are the resulting product
	return _mm256_castpd256_pd128(x3);
}
#else
inline interval_number interval_number::operator*(const interval_number& b) const
{
	// <a0,a1> * <b0,b1>
	__m128d ssg;
	__m128d llhh, lhhl, ip;

	switch ((_mm_movemask_pd(interval) << 2) | _mm_movemask_pd(b.interval))
	{
	case 0: // -+ * -+: <min(<a0*b1>,<a1*b0>), max(<a0*b0>,<a1*b1>)>
		llhh = _mm_mul_pd(interval, b.interval);
		lhhl = _mm_mul_pd(interval, _mm_shuffle_pd(b.interval, b.interval, 1));
		return interval_number(_mm_max_pd(_mm_unpacklo_pd(llhh, lhhl), _mm_unpackhi_pd(llhh, lhhl)));
	case 1: // -+ * --: <b0*a1, b0*a0>
		return interval_number(_mm_mul_pd(_mm_shuffle_pd(b.interval, b.interval, 3), _mm_shuffle_pd(interval, interval, 1)));
	case 2: // -+ * ++: <b1*a0, b1*a1>
		return interval_number(_mm_mul_pd(_mm_shuffle_pd(b.interval, b.interval, 0), interval));
	case 4: // -- * -+: <a0*b1, a0*b0>
		return interval_number(_mm_mul_pd(_mm_shuffle_pd(interval, interval, 3), _mm_shuffle_pd(b.interval, b.interval, 1)));
	case 5: // -- * --: <a1*b1, a0*b0>
		ip = _mm_mul_pd(_mm_xor_pd(interval, sign_high_mask()), b.interval);
		return interval_number(_mm_shuffle_pd(ip, ip, 1));
	case 6: // -- * ++: <a0*b1, a1*b0>
		ssg = _mm_xor_pd(b.interval, sign_low_mask());
		return interval_number(_mm_mul_pd(interval, _mm_shuffle_pd(ssg, ssg, 1)));
	case 8: // ++ * -+: <a1*b0, a1*b1>
		return interval_number(_mm_mul_pd(_mm_shuffle_pd(interval, interval, 0), b.interval));
	case 9: // ++ * --: <b0*a1, b1*a0>
		ssg = _mm_xor_pd(interval, sign_low_mask());
		return interval_number(_mm_mul_pd(b.interval, _mm_shuffle_pd(ssg, ssg, 1)));
	case 10: // ++ * ++: <a0*b0, a1*b1>
		return interval_number(_mm_mul_pd(interval, _mm_xor_pd(b.interval, sign_low_mask())));
	}

	return interval_number(NAN);
}
#endif

inline interval_number interval_number::operator*(const double b) const {
	if (b >= 0) return interval_number(_mm_mul_pd(interval, _mm_set1_pd(b)));
	else return interval_number(_mm_mul_pd(_mm_shuffle_pd(interval, interval, 1), _mm_set1_pd(-b)));
}

inline interval_number interval_number::operator/(const double b) const {
	if (b >= 0) return interval_number(_mm_div_pd(interval, _mm_set1_pd(b)));
	else return interval_number(_mm_div_pd(_mm_shuffle_pd(interval, interval, 1), _mm_set1_pd(-b)));
}

inline interval_number interval_number::abs() const {
	switch (_mm_movemask_pd(interval))
	{
	case 0: // Hi>0, Lo<0
		return _mm_and_pd(_mm_max_pd(interval, _mm_shuffle_pd(interval, interval, 1)), all_high_mask());
	case 1: // Hi<0, Lo<0
		return _mm_shuffle_pd(interval, interval, 1);
	}
	return *this; // If Hi>0, Lo>0 OR invalid interval == case 3 above
}

inline interval_number interval_number::sqr() const {
	const interval_number av = abs();
	return interval_number(_mm_mul_pd(av.getLowSwitched(), av.interval));
}

inline interval_number min(const interval_number& a, const interval_number& b) {
	const __m128d ai = a.getLowSwitched(), bi = b.getLowSwitched();
	const __m128d m = _mm_min_pd(ai, bi);
	return interval_number(m).getLowSwitched();
}

inline interval_number max(const interval_number& a, const interval_number& b) {
	const __m128d ai = a.getLowSwitched(), bi = b.getLowSwitched();
	const __m128d m = _mm_max_pd(ai, bi);
	return interval_number(m).getLowSwitched();
}

inline interval_number interval_number::inverse() const
{
	const int m = _mm_movemask_pd(interval);
	if (m == 1 || m == 2)
	{
		const __m128d den = _mm_shuffle_pd(interval, interval, 1);
		const __m128d frac = _mm_div_pd(minus_one(), den);
		return interval_number(frac);
	}
	else
	{
		return interval_number(NAN);
	}
}

inline interval_number interval_number::pow(unsigned int e) const {
	if (e == 0) return interval_number(1.0);

	__m128d uns = _mm_and_pd(interval, sign_fabs_mask());
	__m128d ui = interval;

	if (!(e & (unsigned int)1)) { // e is even
		const int s = _mm_movemask_pd(interval);
		if (s == 0) {
			__m128d swapped = _mm_shuffle_pd(uns, uns, 1);
			if (_mm_comigt_sd(swapped, uns)) {
				ui = _mm_shuffle_pd(uns, uns, 1);
				uns = swapped;
			}
			uns = _mm_and_pd(uns, all_high_mask());
		}
		else if (s == 1) {
			ui = _mm_shuffle_pd(ui, ui, 1);
			uns = _mm_shuffle_pd(uns, uns, 1);
		}
	}

	while (--e) ui = _mm_mul_pd(ui, uns);
	return interval_number(ui);
}

#else // USE_SIMD_INSTRUCTIONS

inline const double* interval_number::getInterval() const { return (const double*)&min_low; }

inline interval_number::interval_number() { }
inline interval_number::interval_number(double a) : min_low(-a), high(a) {}
inline interval_number::interval_number(double minf, double sup) : min_low(minf), high(sup) {}
inline interval_number::interval_number(const interval_number& b) : min_low(b.min_low), high(b.high) {}

inline double interval_number::minus_inf() const { return min_low; }
inline double interval_number::inf() const { return -min_low; }
inline double interval_number::sup() const { return high; }

inline bool interval_number::isNegative() const { return (high < 0); }
inline bool interval_number::isPositive() const { return (min_low < 0); }
inline void interval_number::negate() { std::swap(min_low, high); }

inline bool interval_number::operator<(const double b) const { return (high < b); }

inline interval_number& interval_number::operator=(const interval_number& b) { min_low = b.min_low; high = b.high; return *this; }

inline interval_number interval_number::operator+(const interval_number& b) const { return interval_number(min_low + b.min_low, high + b.high); }

inline interval_number interval_number::operator-(const interval_number& b) const { return interval_number(b.high + min_low, high + b.min_low); }

inline interval_number interval_number::operator-() const { return interval_number(high, min_low); }
inline interval_number interval_number::operator+(const double b) const { return interval_number(min_low - b, high + b); }
inline interval_number interval_number::operator-(const double b) const { return interval_number(min_low + b, high - b); }

inline interval_number& interval_number::operator+=(const interval_number& b) { min_low += b.min_low; high += b.high; return *this; }
inline interval_number& interval_number::operator-=(const interval_number& b) { return operator=(*this - b); }
inline interval_number& interval_number::operator*=(const interval_number& b) { return operator=(*this * b); }
inline interval_number& interval_number::operator+=(const double b) { min_low -= b; high += b; return *this; }
inline interval_number& interval_number::operator-=(const double b) { min_low += b; high -= b; return *this; }
inline interval_number& interval_number::operator*=(const double b) { return operator=(*this * b); }
inline interval_number& interval_number::operator/=(const double b) { return operator=(*this / b); }

inline interval_number min(const interval_number& a, const interval_number& b) {
	return interval_number(std::max(a.min_low, b.min_low), std::min(a.high, b.high));
}

inline interval_number max(const interval_number& a, const interval_number& b) {
	return interval_number(std::min(a.min_low, b.min_low), std::max(a.high, b.high));
}

inline bool interval_number::operator>(const double b) const { return (min_low < -b); }
inline bool interval_number::operator==(const double b) const { return (high == b && min_low == -b); }

inline int interval_number::sign() const { return (isNegative()) ? (-1) : (1); } // Zero is not accounted for

inline interval_number interval_number::operator*(const interval_number& b) const
{
	typedef union error_approx_type_t
	{
		double d;
		uint64_t u;

		inline error_approx_type_t() {}
		inline error_approx_type_t(double a) : d(a) {}
		inline uint64_t is_negative() const { return u >> 63; }
	} casted_double;

	casted_double l1(min_low), h1(high), l2(b.min_low), h2(b.high);
	uint64_t cfg = (l1.is_negative() << 3) + (h1.is_negative() << 2) + (l2.is_negative() << 1) + (h2.is_negative());

	switch (cfg)
	{
	case 10: return interval_number(min_low * (-b.min_low), high * b.high);
	case 8: return interval_number(high * b.min_low, high * b.high);
	case 9: return interval_number(high * b.min_low, (-min_low) * b.high);
	case 2: return interval_number(min_low * b.high, high * b.high);
	case 0:
		double ll, lh, hl, hh;
		ll = min_low * b.min_low; lh = (min_low * b.high); hl = (high * b.min_low); hh = high * b.high;
		if (hl > lh) lh = hl;
		if (ll > hh) hh = ll;
		return interval_number(lh, hh);
	case 1: return interval_number(high * b.min_low, min_low * b.min_low);
	case 6: return interval_number(min_low * b.high, high * (-b.min_low));
	case 4: return interval_number(min_low * b.high, min_low * b.min_low);
	case 5: return interval_number((-high) * b.high, min_low * b.min_low);
	};

	return interval_number(NAN);
}

inline interval_number interval_number::operator*(const double b) const
{
	if (b >= 0) return interval_number(min_low * b, high * b);
	else return interval_number(high * (-b), min_low * (-b));
}

inline interval_number interval_number::operator/(const double b) const
{
	if (b >= 0) return interval_number(min_low / b, high / b);
	else return interval_number(high / (-b), min_low / (-b));
}

inline interval_number interval_number::abs() const {
	if (min_low < 0) return *this;
	if (high < 0) return interval_number(high, min_low);
	return interval_number(0, std::max(high, min_low));
}

inline interval_number interval_number::sqr() const {
	if (min_low < 0) return interval_number(-min_low * min_low, high * high);
	if (high < 0) return interval_number(-high * high, min_low * min_low);
	if (min_low < high) return interval_number(0, high * high);
	return interval_number(0, min_low * min_low);
}

inline interval_number interval_number::inverse() const
{
	if ((min_low < 0 && high>0) || (min_low > 0 && high < 0))
	{
		return interval_number(-1.0 / high, -1.0 / min_low);
	}
	else
	{
		return interval_number(NAN);
	}
}

inline interval_number interval_number::pow(unsigned int e) const {

	const double _uinf = fabs(min_low), _usup = fabs(high);

	if (e & (unsigned int)1) { // If e is odd
		double _uml = min_low, _uh = high;
		while (--e) { _uml *= _uinf; _uh *= _usup; }
		return interval_number(_uml, _uh);
	}
	else { // e is even
		if (e == 0) return interval_number(1.0);

		if (_uinf > _usup) {
			double _uml = (min_low > 0 && high > 0) ? 0 : (-_usup);
			double _uh = _uinf;
			while (--e) { _uml *= _usup; _uh *= _uinf; }
			return interval_number(_uml, _uh);
		}
		else {
			double _uml = (min_low > 0 && high > 0) ? 0 : (-_uinf);
			double _uh = _usup;
			while (--e) { _uml *= _uinf; _uh *= _usup; }
			return interval_number(_uml, _uh);
		}
	}
}

#endif // USE_SIMD_INSTRUCTIONS

inline double interval_number::width() const { return sup() - inf(); }

inline bool interval_number::signIsReliable() const { return (isNegative() || isPositive()); } // Zero is not accounted for
inline bool interval_number::containsZero() const { return !signIsReliable(); }

inline bool interval_number::isNAN() const { return sup() != sup(); }

inline double interval_number::getMid() const { return (inf() + sup()) / 2; }
inline bool interval_number::isExact() const { return inf() == sup(); }

inline bool interval_number::operator<(const interval_number& b) const { return (sup() < b.inf()); }
inline bool interval_number::operator>(const interval_number& b) const { return (inf() > b.sup()); }

inline bool interval_number::operator<=(const interval_number& b) const { return (sup() <= b.inf()); }
inline bool interval_number::operator>=(const interval_number& b) const { return (inf() >= b.sup()); }

inline bool interval_number::operator==(const interval_number& b) const { return (sup() == inf() && sup() == b.inf() && sup() == b.sup()); }

inline bool interval_number::operator!=(const interval_number& b) const { return operator<(b) || operator>(b); }

inline interval_number sqrt(const interval_number& p)
{
	const double inf = p.inf();
	const double sup = p.sup();
	if (inf < 0 || sup < 0) return interval_number(NAN);
	const double srinf = sqrt(inf);
	const double srsup = sqrt(sup);
	if (srinf * srinf > inf) return interval_number((-nextafter(srinf, 0)), srsup);
	else return interval_number(-srinf, srsup);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// FUNCTION IMPLEMENTATIONS - EXPANSIONS
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

// Allocate extra-memory
//#define AllocDoubles(n) ((double *)malloc((n) * sizeof(double)))
//#define FreeDoubles(p) (free(p))
#define AllocDoubles(n) ((double *)expansionObject::mempool.alloc((n) * sizeof(double)))
#define FreeDoubles(p) (expansionObject::mempool.release(p))

inline void expansionObject::Quick_Two_Sum(const double a, const double b, double& x, double& y) { x = a + b; y = b - (x - a); }

inline void expansionObject::quick_Two_Sum(const double a, const double b, double* xy) { Quick_Two_Sum(a, b, xy[1], xy[0]); }

inline void expansionObject::Two_Sum(const double a, const double b, double& x, double& y) {
	double _bv;
	x = a + b; _bv = x - a; y = (a - (x - _bv)) + (b - _bv);
}

inline void expansionObject::two_Sum(const double a, const double b, double* xy) { Two_Sum(a, b, xy[1], xy[0]); }

inline void expansionObject::Two_One_Sum(const double a1, const double a0, const double b, double& x2, double& x1, double& x0) {
	double _i;
	Two_Sum(a0, b, _i, x0); Two_Sum(a1, _i, x2, x1);
}

inline void expansionObject::two_One_Sum(const double *a, const double b, double *x) {
	Two_One_Sum(a[1], a[0], b, x[2], x[1], x[0]);
}

inline void expansionObject::Two_Diff(const double a, const double b, double& x, double& y) {
	double _bv;
	x = a - b; _bv = a - x; y = (a - (x + _bv)) + (_bv - b);
}

inline void expansionObject::Two_One_Diff(const double a1, const double a0, const double b, double& x2, double& x1, double& x0) {
	double _i;
	Two_Diff(a0, b, _i, x0); Two_Sum(a1, _i, x2, x1);
}

inline void expansionObject::two_Diff(const double a, const double b, double* xy) { Two_Diff(a, b, xy[1], xy[0]); }

// Products
#ifndef USE_AVX2_INSTRUCTIONS 
inline void expansionObject::Split(double a, double& _ah, double& _al) {
	double _c = 1.3421772800000003e+008 * a;
	_ah = _c - (_c - a); _al = a - _ah;
}

inline void expansionObject::Two_Prod_PreSplit(double a, double b, double _bh, double _bl, double& x, double& y) {
	double _ah, _al;
	x = a * b;
	Split(a, _ah, _al);
	y = (_al * _bl) - (((x - (_ah * _bh)) - (_al * _bh)) - (_ah * _bl));
}

inline void expansionObject::Two_Product_2Presplit(double a, double _ah, double _al, double b, double _bh, double _bl, double& x, double& y) {
	x = a * b; y = (_al * _bl) - (((x - _ah * _bh) - (_al * _bh)) - (_ah * _bl));
}
#endif

// [x,y] = [a]*[b]		 Multiplies two expansions [a] and [b] of length one
inline void expansionObject::Two_Prod(const double a, const double b, double* xy) { Two_Prod(a, b, xy[1], xy[0]); }


// [x,y] = [a]^2		Squares an expansion of length one
inline void expansionObject::Square(const double a, double* xy) { Square(a, xy[1], xy[0]); }

// [x2,x1,x0] = [a1,a0]-[b]		Subtracts an expansion [b] of length one from an expansion [a1,a0] of length two
inline void expansionObject::two_One_Diff(const double a1, const double a0, const double b, double& x2, double& x1, double& x0)
{
	Two_One_Diff(a1, a0, b, x2, x1, x0);
}

inline void expansionObject::two_One_Diff(const double* a, const double b, double* x) { two_One_Diff(a[1], a[0], b, x[2], x[1], x[0]); }

// [x3,x2,x1,x0] = [a1,a0]*[b]		Multiplies an expansion [a1,a0] of length two by an expansion [b] of length one
inline void expansionObject::Two_One_Prod(const double* a, const double b, double* x) { Two_One_Prod(a[1], a[0], b, x[3], x[2], x[1], x[0]); }

// [x3,x2,x1,x0] = [a1,a0]+[b1,b0]		Calculates the sum of two expansions of length two
inline void expansionObject::Two_Two_Sum(const double a1, const double a0, const double b1, const double b0, double& x3, double& x2, double& x1, double& x0) {
	double _j, _0;
	Two_One_Sum(a1, a0, b0, _j, _0, x0); Two_One_Sum(_j, _0, b1, x3, x2, x1);
}

inline void expansionObject::Two_Two_Sum(const double* a, const double* b, double* xy) { Two_Two_Sum(a[1], a[0], b[1], b[0], xy[3], xy[2], xy[1], xy[0]); }

// [x3,x2,x1,x0] = [a1,a0]-[b1,b0]		Calculates the difference between two expansions of length two
inline void expansionObject::Two_Two_Diff(const double a1, const double a0, const double b1, const double b0, double& x3, double& x2, double& x1, double& x0) {
	double _j, _0, _u3;
	Two_One_Diff(a1, a0, b0, _j, _0, x0); Two_One_Diff(_j, _0, b1, _u3, x2, x1); x3 = _u3;
}

inline void expansionObject::Two_Two_Diff(const double* a, const double* b, double* x) { Two_Two_Diff(a[1], a[0], b[1], b[0], x[3], x[2], x[1], x[0]); }

// Calculates the second component 'y' of the expansion [x,y] = [a]-[b] when 'x' is known
inline void expansionObject::Two_Diff_Back(const double a, const double b, double& x, double& y) {
	double _bv;
	_bv = a - x; y = (a - (x + _bv)) + (_bv - b);
}

inline void expansionObject::Two_Diff_Back(const double a, const double b, double* xy) { Two_Diff_Back(a, b, xy[1], xy[0]); }

inline void expansionObject::Two_Two_Prod(const double* a, const double* b, double* xy) { Two_Two_Prod(a[1], a[0], b[1], b[0], xy); }

// [e] = -[e]		Inplace inversion
inline void expansionObject::Gen_Invert(const int elen, double* e) { for (int i = 0; i < elen; i++) e[i] = -e[i]; }

// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
inline int expansionObject::Gen_Sum_With_Alloc(const int elen, const double* e, const int flen, const double* f, double** h)
{
	*h = AllocDoubles(elen + flen);
	return Gen_Sum(elen, e, flen, f, *h);
}

// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
inline int expansionObject::Gen_Diff_With_Alloc(const int elen, const double* e, const int flen, const double* f, double** h)
{
	*h = AllocDoubles(elen + flen);
	return Gen_Diff(elen, e, flen, f, *h);
}

// [h] = [e] * 2		Multiplies an expansion by 2
// 'h' must be allocated by the caller with at least elen components. This is exact up to overflows.
inline void expansionObject::Double(const int elen, const double* e, double* h) { for (int i = 0; i < elen; i++) h[i] = 2 * e[i]; }

// [h] = [e] * n		Multiplies an expansion by n
// If 'n' is a power of two, the multiplication is exact
inline void expansionObject::ExactScale(const int elen, double* e, const double n) { for (int i = 0; i < elen; i++) e[i] *= n; }

inline void expansionObject::print(const int elen, const double* e) { for (int i = 0; i < elen; i++) printf("%e ", e[i]); printf("\n"); }

#ifdef USE_AVX2_INSTRUCTIONS
inline void vfast_Two_Sum(__m128d a, __m128d b, __m128d& x, __m128d& y) {
	x = _mm_add_sd(a, b);
	y = _mm_sub_sd(b, _mm_sub_sd(x, a));
}

inline void vtwo_Sum(__m128d a, __m128d b, __m128d& x, __m128d& y) {
	x = _mm_add_sd(a, b);
	__m128d vbv = _mm_sub_sd(x, a);
	y = _mm_add_sd(_mm_sub_sd(a, _mm_sub_sd(x, vbv)), _mm_sub_sd(b, vbv));
}

inline void vtwo_Prod(__m128d a, __m128d b, __m128d& x, __m128d& y) {
	x = _mm_mul_sd(a, b);
	y = _mm_fmsub_sd(a, b, x);
}
#endif

inline void expansionObject::Two_Prod(const double a, const double b, double& x, double& y)
{
#ifdef USE_AVX2_INSTRUCTIONS
	__m128d x1, x2, x3;
	x1 = _mm_set_sd(a);
	x2 = _mm_set_sd(b);
	vtwo_Prod(x1, x2, x3, x2);
	x = _mm_cvtsd_f64(x3);
	y = _mm_cvtsd_f64(x2);
#else
	double _ah, _al, _bh, _bl;
	x = a * b;
	Split(a, _ah, _al); Split(b, _bh, _bl);
	y = ((_ah * _bh - x) + _ah * _bl + _al * _bh) + _al * _bl;
#endif
}

inline void expansionObject::Square(const double a, double& x, double& y)
{
#ifdef USE_AVX2_INSTRUCTIONS
	__m128d x1, x2, x3;
	x1 = _mm_set_sd(a);
	vtwo_Prod(x1, x1, x3, x2);
	x = _mm_cvtsd_f64(x3);
	y = _mm_cvtsd_f64(x2);
#else
	double _ah, _al;
	x = a * a;
	Split(a, _ah, _al);
	y = (_al * _al) - ((x - (_ah * _ah)) - ((_ah + _ah) * _al));
#endif
}

inline void expansionObject::Two_One_Prod(const double a1, const double a0, const double b, double& x3, double& x2, double& x1, double& x0)
{
#ifdef USE_AVX2_INSTRUCTIONS
	double _i, _j, _k, _0;
	Two_Prod(a0, b, _i, x0); Two_Prod(a1, b, _j, _0);
	Two_Sum(_i, _0, _k, x1); Quick_Two_Sum(_j, _k, x3, x2);
#else
	double _bh, _bl, _i, _j, _0, _k;
	Split(b, _bh, _bl);
	Two_Prod_PreSplit(a0, b, _bh, _bl, _i, x0); Two_Prod_PreSplit(a1, b, _bh, _bl, _j, _0);
	Two_Sum(_i, _0, _k, x1); Quick_Two_Sum(_j, _k, x3, x2);
#endif
}

inline void expansionObject::Two_Two_Prod(const double a1, const double a0, const double b1, const double b0, double* h)
{
#ifdef USE_AVX2_INSTRUCTIONS
	double _m, _n, _i, _j, _k, _0, _1, _2, _l;
	Two_Prod(a0, b0, _i, h[0]);
	Two_Prod(a1, b0, _j, _0);
	Two_Sum(_i, _0, _k, _1);
	Quick_Two_Sum(_j, _k, _l, _2);
	Two_Prod(a0, b1, _i, _0);
	Two_Sum(_1, _0, _k, h[1]);
	Two_Sum(_2, _k, _j, _1);
	Two_Sum(_l, _j, _m, _2);
	Two_Prod(a1, b1, _j, _0);
	Two_Sum(_i, _0, _n, _0);
	Two_Sum(_1, _0, _i, h[2]);
	Two_Sum(_2, _i, _k, _1);
	Two_Sum(_m, _k, _l, _2);
	Two_Sum(_j, _n, _k, _0);
	Two_Sum(_1, _0, _j, h[3]);
	Two_Sum(_2, _j, _i, _1);
	Two_Sum(_l, _i, _m, _2);
	Two_Sum(_1, _k, _i, h[4]);
	Two_Sum(_2, _i, _k, h[5]);
	Two_Sum(_m, _k, h[7], h[6]);
#else
	double _ah, _al, _bh, _bl, _ch, _cl, _m, _n, _i, _j, _0, _1, _2, _k, _l;
	Split(a0, _ah, _al);
	Split(b0, _bh, _bl);
	Two_Product_2Presplit(a0, _ah, _al, b0, _bh, _bl, _i, h[0]);
	Split(a1, _ch, _cl);
	Two_Product_2Presplit(a1, _ch, _cl, b0, _bh, _bl, _j, _0);
	Two_Sum(_i, _0, _k, _1);
	Quick_Two_Sum(_j, _k, _l, _2);
	Split(b1, _bh, _bl);
	Two_Product_2Presplit(a0, _ah, _al, b1, _bh, _bl, _i, _0);
	Two_Sum(_1, _0, _k, h[1]);
	Two_Sum(_2, _k, _j, _1);
	Two_Sum(_l, _j, _m, _2);
	Two_Product_2Presplit(a1, _ch, _cl, b1, _bh, _bl, _j, _0);
	Two_Sum(_i, _0, _n, _0);
	Two_Sum(_1, _0, _i, h[2]);
	Two_Sum(_2, _i, _k, _1);
	Two_Sum(_m, _k, _l, _2);
	Two_Sum(_j, _n, _k, _0);
	Two_Sum(_1, _0, _j, h[3]);
	Two_Sum(_2, _j, _i, _1);
	Two_Sum(_l, _i, _m, _2);
	Two_Sum(_1, _k, _i, h[4]);
	Two_Sum(_2, _i, _k, h[5]);
	Two_Sum(_m, _k, h[7], h[6]);
#endif
}

inline int expansionObject::Gen_Sum(const int elen, const double* e, const int flen, const double* f, double* h)
{
	double Q, Qn, hh, s;
	const double* en = e, * fn = f, * elast = e + elen, * flast = f + flen;
	int h_k = 0;

	Q = (fabs(*fn) > fabs(*en)) ? (*en++) : (*fn++);

	if ((en < elast) && (fn < flast))
	{
		s = (fabs(*fn) > fabs(*en)) ? (*en++) : (*fn++);
		Quick_Two_Sum(s, Q, Qn, hh);
		Q = Qn;
		if (hh != 0.0) h[h_k++] = hh;
		while ((en < elast) && (fn < flast))
		{
			s = (fabs(*fn) > fabs(*en)) ? (*en++) : (*fn++);
			Two_Sum(s, Q, Qn, hh);
			Q = Qn;
			if (hh != 0.0) h[h_k++] = hh;
		}
	}

	while (en < elast)
	{
		Two_Sum(Q, (*en), Qn, hh);
		Q = Qn;
		if (hh != 0.0) h[h_k++] = hh;
		en++;
	}

	while (fn < flast)
	{
		Two_Sum(Q, (*fn), Qn, hh);
		Q = Qn;
		if (hh != 0.0) h[h_k++] = hh;
		fn++;
	}
	if ((Q != 0.0) || (h_k == 0)) h[h_k++] = Q;

	return h_k;
}

inline int expansionObject::Gen_Diff(const int elen, const double* e, const int flen, const double* f, double* h)
{
	double Q, Qn, hh, s;
	const double* en = e, * fn = f, * elast = e + elen, * flast = f + flen;
	int h_k = 0;

	Q = (fabs(*fn) > fabs(*en)) ? (*en++) : (-*fn++);

	if ((en < elast) && (fn < flast))
	{
		s = (fabs(*fn) > fabs(*en)) ? (*en++) : (-*fn++);
		Quick_Two_Sum(s, Q, Qn, hh);
		Q = Qn;
		if (hh != 0.0) h[h_k++] = hh;
		while ((en < elast) && (fn < flast))
		{
			s = (fabs(*fn) > fabs(*en)) ? (*en++) : (-*fn++);
			Two_Sum(s, Q, Qn, hh);
			Q = Qn;
			if (hh != 0.0) h[h_k++] = hh;
		}
	}

	while (en < elast)
	{
		Two_Sum(Q, (*en), Qn, hh);
		Q = Qn;
		if (hh != 0.0) h[h_k++] = hh;
		en++;
	}

	while (fn < flast)
	{
		s = *fn++;
		Two_Diff(Q, s, Qn, hh);
		Q = Qn;
		if (hh != 0.0) h[h_k++] = hh;
	}
	if ((Q != 0.0) || (h_k == 0)) h[h_k++] = Q;

	return h_k;
}


inline int expansionObject::Gen_Scale(const int elen, const double* e, const double b, double* h)
{
	double Q, sum, hh, pr1, pr0;
	const double* ei = e, * elast = e + elen;

	int k = 0;
	Two_Prod(*ei, b, Q, hh);
	if (hh != 0) h[k++] = hh;

	while (++ei < elast) {
		Two_Prod(*ei, b, pr1, pr0);
		Two_Sum(Q, pr0, sum, hh);
		if (hh != 0) h[k++] = hh;
		Quick_Two_Sum(pr1, sum, Q, hh);
		if (hh != 0) h[k++] = hh;
	}
	if ((Q != 0.0) || (k == 0)) h[k++] = Q;
	return k;
}


inline void expansionObject::Two_Square(const double& a1, const double& a0, double* x)
{
	double _j, _0, _k, _1, _l, _2;
	Square(a0, _j, x[0]);
	_0 = a0 + a0;
	Two_Prod(a1, _0, _k, _1);
	Two_One_Sum(_k, _1, _j, _l, _2, x[1]);
	Square(a1, _j, _1);
	Two_Two_Sum(_j, _1, _l, _2, x[5], x[4], x[3], x[2]);
}

inline void expansionObject::Two_Square(const double* a, double* x) { Two_Square(a[1], a[0], x); }

inline int expansionObject::Sub_product(const int alen, const double* a, const int blen, const double* b, double* h)
{
	if (alen == 1) return Gen_Scale(blen, b, a[0], h);
	int partial = 2 * alen * blen;
	int allmem = 2 * (partial + blen);
	double ph1_p[1024];
	double* ph1 = (allmem > 1024) ? (AllocDoubles(allmem)) : (ph1_p);
	double* ph2 = ph1 + partial;
	double* th = ph2 + partial;
	double* ph[2] = { ph1, ph2 };
	int first = 0;
	int phl = Gen_Scale(blen, b, a[0], ph[0]);

	for (int i = 1; i < alen; i++)
	{
		int thl = Gen_Scale(blen, b, a[i], th);
		first = i & 1;
		phl = Gen_Sum(phl, ph[(i + 1) & 1], thl, th, ph[first]);
	}
	if (first) for (int i = 0; i < phl; i++) h[i] = ph2[i];
	else for (int i = 0; i < phl; i++) h[i] = ph1[i];
	if (allmem > 1024) FreeDoubles(ph1);
	return phl;
}


inline int expansionObject::Gen_Product(const int alen, const double* a, const int blen, const double* b, double* h)
{
	if (blen == 1) return Gen_Scale(alen, a, b[0], h);
	else if (alen < blen) return Sub_product(alen, a, blen, b, h);
	else return Sub_product(blen, b, alen, a, h);
}


inline double expansionObject::To_Double(const int elen, const double* e)
{
	double Q = e[0];
	for (int e_i = 1; e_i < elen; e_i++) Q += e[e_i];
	return Q;
}

inline int expansionObject::Gen_Product_With_Alloc(const int alen, const double* a, const int blen, const double* b, double** h)
{
	int h_len = alen * blen * 2;
	if (h_len < 8) h_len = 8;
	*h = AllocDoubles(h_len);
	return Gen_Product(alen, a, blen, b, *h);
}

inline int expansionObject::Double_With_PreAlloc(const int elen, const double* e, double** h, const int hlen)
{
	int newlen = elen;
	if (hlen < newlen) *h = AllocDoubles(newlen);
	Double(elen, e, *h);
	return newlen;
}

inline int expansionObject::Gen_Scale_With_PreAlloc(const int elen, const double* e, const double& b, double** h, const int hlen)
{
	int newlen = elen * 2;
	if (hlen < newlen) *h = AllocDoubles(newlen);
	return Gen_Scale(elen, e, b, *h);
}

inline int expansionObject::Gen_Sum_With_PreAlloc(const int elen, const double* e, const int flen, const double* f, double** h, const int hlen)
{
	int newlen = elen + flen;
	if (hlen < newlen) *h = AllocDoubles(newlen);
	return Gen_Sum(elen, e, flen, f, *h);
}

inline int expansionObject::Gen_Diff_With_PreAlloc(const int elen, const double* e, const int flen, const double* f, double** h, const int hlen)
{
	int newlen = elen + flen;
	if (hlen < newlen) *h = AllocDoubles(newlen);
	return Gen_Diff(elen, e, flen, f, *h);
}

inline int expansionObject::Gen_Product_With_PreAlloc(const int alen, const double* a, const int blen, const double* b, double** h, const int hlen)
{
	int newlen = alen * blen * 2;
	if (hlen < newlen || hlen < 8)
	{
		if (newlen < 8) newlen = 8;
		*h = AllocDoubles(newlen);
	}
	return Gen_Product(alen, a, blen, b, *h);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// FUNCTION IMPLEMENTATIONS - BIGNATURAL
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef USE_GNU_GMP_CLASSES

inline uint32_t* bignatural::BN_ALLOC(uint32_t num_bytes) { return (uint32_t*)nfgMemoryPool.alloc(num_bytes); }
inline void bignatural::BN_FREE(uint32_t* ptr) { nfgMemoryPool.release(ptr); }
inline bignatural::bignatural() : m_capacity(0), m_size(0), digits(NULL) { }
inline bignatural::~bignatural() { BN_FREE(digits); }
inline bignatural::bignatural(const bignatural& m) { init(m); }
inline bignatural::bignatural(bignatural&& m) noexcept : m_capacity(m.m_capacity), m_size(m.m_size), digits(m.digits) {
	m.digits = nullptr;
	m.m_size = m.m_capacity = 0;
}

inline bignatural::bignatural(uint32_t m) { init(m); }

inline bignatural::bignatural(uint64_t m) { init(m); }

inline bignatural::bignatural(FILE* f) { init(f); }

inline const uint32_t& bignatural::back() const { return digits[m_size - 1]; }

inline const uint32_t& bignatural::operator[](int i) const { return digits[i]; }

inline uint32_t bignatural::size() const { return m_size; }

inline bool bignatural::empty() const { return m_size == 0; }

inline bignatural bignatural::operator<<(uint32_t n) const { bignatural t(*this); t <<= n; return t; }

inline bignatural bignatural::operator>>(uint32_t n) const { bignatural t(*this); t >>= n; return t; }

inline bool bignatural::operator<=(const bignatural& b) const { return b >= *this; }
inline bool bignatural::operator<(const bignatural& b) const { return b > *this; }

inline bignatural bignatural::operator+(const bignatural& b) const { bignatural n(*this); return (n += b); };
inline bignatural bignatural::operator+(const uint32_t b) const { bignatural n(*this); return (n += b); }
inline bignatural bignatural::operator+(const uint64_t b) const { return operator+(bignatural(b)); }

inline bignatural bignatural::operator-(const bignatural& b) const { bignatural n(*this); return (n -= b); };

inline bignatural& bignatural::operator*=(const bignatural& b) { operator=(*this * b); return *this; }

inline void bignatural::operator|=(uint32_t i) { if (m_size) digits[m_size - 1] |= i; else operator=(i); }

inline uint32_t& bignatural::back() { return digits[m_size - 1]; }

inline void bignatural::pop_back() { m_size--; }

inline uint32_t& bignatural::operator[](int i) { return digits[i]; }

inline void bignatural::reserve(uint32_t n) { if (n > m_capacity) increaseCapacity(n); }

inline void bignatural::resize(uint32_t n) { reserve(n); m_size = n; }

inline void bignatural::fill(uint32_t v) {
	memset(digits, (int)v, m_size << 2);
}

inline uint32_t bignatural::countEndingZeroesLSL() const {
	return (uint32_t)nfg_count_rz(back());
}

inline size_t bignatural::scan_uint64_t(FILE* fp, uint64_t& t) {
	const uint64_t limit_n = (UINT64_MAX - 9) / 10;
	uint64_t n = 0;
	size_t i = 0;

	// Quick scan of the first two limbs
	int c = fgetc(fp);
	while (isdigit(c) && n <= limit_n) {
		i++;
		n *= 10;
		n += (c - '0');
		c = fgetc(fp);
	}
	if (c != EOF) ungetc(c, fp);

	t = n;
	return i;
}

inline void bignatural::init(FILE* fp) {
	uint64_t n, m;
	scan_uint64_t(fp, n);
	init(n);
	size_t r;
	while ((r = scan_uint64_t(fp, n)) != 0) {
		m = 1;
		while (r--) m *= 10;
		operator*=(m);
		operator+=(n);
	}
}

inline void bignatural::init(const bignatural& m) {
	m_size = m.m_size;
	m_capacity = m.m_capacity;
	if (m_capacity) {
		digits = (uint32_t*)BN_ALLOC(sizeof(uint32_t) * m_capacity);
		if (m_size) memcpy(digits, m.digits, sizeof(uint32_t) * m_size);
	}
	else digits = NULL;
}

inline void bignatural::init(const uint64_t m) {
	if (m == 0) {
		m_size = m_capacity = 0;
		digits = NULL;
	}
	else if (m <= UINT32_MAX) {
		m_size = m_capacity = 1;
		digits = (uint32_t*)BN_ALLOC(sizeof(uint32_t));
		digits[0] = (uint32_t)m;
	}
	else {
		m_size = m_capacity = 2;
		digits = (uint32_t*)BN_ALLOC(sizeof(uint32_t) * 2);
		digits[0] = (uint32_t)(m >> 32);
		digits[1] = (uint32_t)(m);
	}
}

inline void bignatural::init(const uint32_t m) {
	if (m == 0) {
		m_size = m_capacity = 0;
		digits = NULL;
	}
	else {
		m_size = m_capacity = 1;
		digits = (uint32_t*)BN_ALLOC(sizeof(uint32_t));
		digits[0] = m;
	}
}

inline bool bignatural::toUint64(uint64_t& n) const {
	if (m_size == 0) n = 0;
	else if (m_size == 1) n = digits[0];
	else if (m_size == 2) { n = (((uint64_t)digits[0]) << 32) + digits[1]; }
	else return false;

	return true;
}

inline bool bignatural::toUint32(uint32_t& n) const {
	if (m_size == 0) n = 0;
	else if (m_size == 1) n = digits[0];
	else return false;

	return true;
}

inline bignatural& bignatural::operator=(const bignatural& m) {
	if (digits != m.digits) {
		if (m_capacity >= m.m_size) {
			m_size = m.m_size;
			if (m_size) memcpy(digits, m.digits, m_size << 2);
		}
		else {
			BN_FREE(digits);
			init(m);
		}
	}

	return *this;
}

inline bignatural& bignatural::operator=(const uint64_t m) {
	BN_FREE(digits);
	init(m);
	return *this;
}

inline void bignatural::operator<<=(uint32_t n) {
	if (m_size == 0) return;
	leftShift(n);
}

inline void bignatural::leftShift(uint32_t n) {
	uint32_t s = n & 0x0000001f;
	uint32_t lz = countLeadingZeroes();
	uint32_t s2 = 32 - s;

	if (lz < s) { // Need a further limb
		push_back(0);
		auto de = digits + m_size;
		while (digits != --de) {
			*de >>= s2;
			*de |= ((*(de - 1)) << s);
		}
		*de >>= s2;
	}
	else if (s) { // Leading zeroes are enough
		auto dp = digits, de = digits + m_size - 1;
		while (dp != de) {
			*dp <<= s;
			*dp |= ((*(dp + 1)) >> s2);
			dp++;
		}
		*dp <<= s;
	}

	while (n >= 32) {
		push_back(0);
		n -= 32;
	}
}

inline void bignatural::operator>>=(uint32_t n) {
	if ((n >> 6) >= m_size) {
		m_size = 0;
		return;
	}

	while (n >= 32) {
		pop_back();
		n -= 32;
	}
	if (!n) return;

	auto dp = digits, de = digits + m_size;
	while (dp != --de) {
		*de >>= n;
		*de |= ((*(de - 1)) << (32 - n));
	}
	if ((*de >>= n) == 0) pop_front();
}

inline bool bignatural::operator==(const bignatural& b) const {
	if (size() != b.size()) return false;
	auto dp = digits, de = digits + m_size, db = b.digits;
	while (dp != de && *dp == *db) { dp++; db++; }
	return (dp == de);
}

inline bool bignatural::operator!=(const bignatural& b) const {
	if (size() != b.size()) return true;
	auto dp = digits, de = digits + m_size, db = b.digits;
	while (dp != de && *dp == *db) { dp++; db++; }
	return (dp != de);
}

inline bool bignatural::operator>=(const bignatural& b) const {
	const int s = (size() > b.size()) - (size() < b.size());
	if (s) return (s > 0);

	auto dp = digits, de = digits + m_size, db = b.digits;
	while (dp != de && *dp == *db) { dp++; db++; }
	return (dp == de || *dp > *db);
}

inline bool bignatural::operator>(const bignatural& b) const {
	const int s = (size() > b.size()) - (size() < b.size());
	if (s) return (s > 0);

	auto dp = digits, de = digits + m_size, db = b.digits;
	while (dp != de && *dp == *db) { dp++; db++; }
	return (dp != de && *dp > *db);
}

inline bignatural bignatural::operator*(const bignatural& b) const {
	bignatural result;
	result.toProd(*this, b);
	return result;
}

inline bignatural bignatural::operator|(const bignatural& b) const {
	bignatural result;

	if (m_size >= b.size()) {
		result = *this;
		uint32_t* rd = result.digits + m_size;
		uint32_t* bd = b.digits + b.m_size;
		do { *--rd |= *--bd; } while (bd != b.digits);
	}
	else {
		result = b;
		uint32_t* rd = result.digits + m_size;
		uint32_t* bd = b.digits + b.m_size;
		do { *--rd |= *--bd; } while (rd != result.digits);
	}

	return result;
}

inline bignatural bignatural::divide_by(const uint32_t D, uint32_t& remainder) const {
	assert(D != 0);
	//if (D == 0) ip_error("Division by zero\n");
	//if (m_size == 0) return 0;

	// If both dividend fits into 64 bits, use hardware division
	uint64_t n;
	if (toUint64(n)) {
		remainder = n % D;
		return n / D;
	}

	bignatural Q;
	uint32_t next_digit = 0;
	uint64_t dividend = digits[next_digit++];
	for (;;) {
		uint64_t tmp_div = dividend / D;
		if (!Q.empty() || tmp_div) Q.push_back((uint32_t)tmp_div);
		dividend -= (tmp_div * D);
		if (next_digit < m_size) {
			dividend <<= 32;
			dividend += digits[next_digit++];
		}
		else break;
	}
	remainder = (uint32_t)dividend;

	return Q;
}

inline uint32_t bignatural::getNumSignificantBits() const {
	if (!m_size) return 0;
	return (m_size * 32) - nfg_count_lz(digits[0]);
}

inline bool bignatural::getBit(uint32_t b) const {
	const uint32_t dig = (m_size - (b >> 5)) - 1;
	const uint32_t bit = b & 31;
	return (digits[dig] & (1 << bit));
}

// Long division
inline bignatural bignatural::divide_by(const bignatural& divisor, bignatural& remainder) const {
	if (divisor.empty()) ip_error("Division by zero\n");
	if (empty()) return (uint32_t)0;

	// If divisor fits into 32 bits, revert to short division
	uint32_t d32, rem;
	if (divisor.toUint32(d32)) {
		bignatural q = divide_by(d32, rem);
		remainder = rem;
		return q;
	}

	// If both dividend and divisor fit into 64 bits, use hardware division
	uint64_t n, d;
	if (toUint64(n) && divisor.toUint64(d)) {
		remainder = n % d;
		return n / d;
	}

	// If divisor is greater than dividend...
	if (divisor > *this) {
		remainder = *this;
		return (uint32_t)0;
	}

	// Use binary (per-bit) long division
	const bignatural& dividend = *this;

	// Possible optimizations:
	// Pre-allocate number of bits in quotient (bits dividend - bits divisor + 1)

	bignatural quotient, loc_dividend;
	uint32_t next_dividend_bit = dividend.getNumSignificantBits();

	do {
		loc_dividend.push_bit_back(dividend.getBit(--next_dividend_bit));
		if (loc_dividend >= divisor) {
			loc_dividend -= divisor;
			quotient.push_bit_back(1);
		}
		else if (!quotient.empty()) quotient.push_bit_back(0);
	} while (next_dividend_bit);

	remainder = loc_dividend;

	return quotient;
}

inline bignatural bignatural::sqrt() const
{
	bignatural low, tmp, high = *this, mid = *this;

	mid >>= 1U;
	mid += 1U;

	while (high > low)
	{
		if (mid * mid > *this) {
			high = mid - 1U;
			tmp = low + 1U;
			mid -= tmp;
			mid >>= 1;
			mid += tmp;
		}
		else {
			low = mid;
			mid = high - low;
			mid >>= 1;
			mid += (low + 1U);
		}
	}

	return low;
}

inline void bignatural::increaseCapacity(uint32_t new_capacity) {
	m_capacity = new_capacity;
	uint32_t* tmp_d = (uint32_t*)BN_ALLOC(sizeof(uint32_t) * m_capacity);
	if (m_size) memcpy(tmp_d, digits, sizeof(uint32_t) * m_size);
	BN_FREE(digits);
	digits = tmp_d;
}

// Add one most significant digit to this number and make room for it if necessary
inline void bignatural::addOneMostSignificantDigit(uint32_t d) {
	if (m_capacity == m_size) {
		m_capacity++;
		uint32_t* tmp_d = (uint32_t*)BN_ALLOC(sizeof(uint32_t) * m_capacity);
		if (m_size) memcpy(tmp_d + 1, digits, sizeof(uint32_t) * m_size);
		BN_FREE(digits);
		digits = tmp_d;
		digits[0] = d;
	}
	else { // m_capacity > m_size
		uint32_t* td = digits + m_size;
		while (--td != digits) *(td + 1) = *td;
		*(td + 1) = *td;
		*td = d;
	}
	m_size++;
}

// Greatest common divisor (Euclidean algorithm)
inline bignatural bignatural::GCD(const bignatural& D) const {
	bignatural A = *this;
	bignatural B = D;
	bignatural R;
	while (!A.empty() && !B.empty()) {
		A.divide_by(B, R);
		A = B;
		B = R;
	}
	if (A.empty()) return B;
	else return A;
}

// String representation in decimal form
inline std::string bignatural::get_dec_str() const {
	std::string st;
	bignatural N = *this;
	uint32_t R;
	if (N.empty()) return "0";
	while (!N.empty()) {
		N = N.divide_by(10, R);
		st += ('0' + (char)R);
	}
	std::reverse(st.begin(), st.end());

	return st;
}

// String representation in binary form
inline std::string bignatural::get_str() const {
	std::string st;
	char s[33];
	s[32] = 0;
	for (uint32_t j = 0; j < m_size; j++) {
		for (int i = 0; i < 32; i++)
			s[i] = (digits[j] & (((uint32_t)1) << (31 - i))) ? '1' : '0';
		st += s;
	}
	return st;
}

// Count number of zeroes on the right (least significant binary digits)
inline uint32_t bignatural::countEndingZeroes() const {
	if (m_size == 0) return 0;
	uint32_t i = m_size;
	uint32_t shft = 0;
	while (!digits[--i])  shft += 32;

	return shft + nfg_count_rz(digits[i]);
}

inline uint32_t bignatural::countLeadingZeroes() const {
	return (uint32_t)nfg_count_lz(digits[0]);
}

inline void bignatural::toSum(const bignatural& a, const bignatural& b) {
	operator=(a);
	operator+=(b);
}

inline bignatural& bignatural::operator+=(const bignatural& b) {
	if (&b == this) { // If b is this same number just multiply by two
		operator<<=(1U);
		return *this;
	}

	if (m_size == 0) return operator=(b);
	else if (b.m_size == 0) return *this;
	else {
		const uint32_t a_s = m_size;
		const uint32_t b_s = b.m_size;
		uint64_t carry = 0;
		uint32_t* dig_b = b.digits + b_s;

		if (a_s > b_s) {
			uint32_t* dig_a = digits + a_s;
			do {
				const uint64_t da = *(--dig_a);
				const uint64_t db = *(--dig_b);
				const uint64_t sm = da + db + carry;
				*(dig_a) = (uint32_t)sm;
				carry = (sm >> 32);
			} while (dig_b != b.digits);
			do {
				const uint64_t da = *(--dig_a);
				const uint64_t sm = da + carry;
				*(dig_a) = (uint32_t)sm;
				carry = (sm >> 32);
			} while (dig_a != digits);
		}
		else if (a_s == b_s) {
			uint32_t* dig_a = digits + a_s;
			do {
				const uint64_t da = *(--dig_a);
				const uint64_t db = *(--dig_b);
				const uint64_t sm = da + db + carry;
				*(dig_a) = (uint32_t)sm;
				carry = (sm >> 32);
			} while (dig_b != b.digits);
		}
		else { // if (a_s < b_s)
			resize(b_s);
			uint32_t* dig_a = digits + a_s;
			uint32_t* dig_r = digits + b_s;
			do {
				const uint64_t da = *(--dig_a);
				const uint64_t db = *(--dig_b);
				const uint64_t sm = da + db + carry;
				*(--dig_r) = (uint32_t)sm;
				carry = (sm >> 32);
			} while (dig_a != digits);
			do {
				const uint64_t db = *(--dig_b);
				const uint64_t sm = db + carry;
				*(--dig_r) = (uint32_t)sm;
				carry = (sm >> 32);
			} while (dig_b != b.digits);
		}

		if (carry) addOneMostSignificantDigit((uint32_t)carry);
	}
	return *this;
}

inline void bignatural::toDiff(const bignatural& a, const bignatural& b) {
	operator=(a);
	operator-=(b);
}

inline bignatural& bignatural::operator-=(const bignatural& b) {
	assert(operator>=(b));

	if (b.m_size) {
		uint64_t da, db, debt = 0;
		uint32_t* dig_a = digits + m_size;
		const uint32_t* dig_b = b.digits + b.m_size;

		do {
			da = *(--dig_a);
			db = *(--dig_b) + debt;
			debt = (da < db);
			da += (debt << 32);
			*(dig_a) = (uint32_t)(da - db);
		} while (dig_b != b.digits);

		while (dig_a != digits) {
			da = *(--dig_a);
			db = debt;
			debt = (da < db);
			da += (debt << 32);
			*(dig_a) = (uint32_t)(da - db);
		}

		pack();
	}

	return *this;
}

inline void bignatural::toProd(const bignatural& a, const bignatural& b) {
	assert(m_size == 0 && m_capacity == 0); // This assumes that the number is zero!
	if (a.empty()) operator=(a);
	else if (b.empty()) operator=(b);
	else {
		// Uses the naive multiplication. This is the best choice for not-too-large factors
		// Consider implementing Karatsuba's algorithm for larger numbers
		m_size = m_capacity = (a.m_size + b.m_size);
		digits = (uint32_t*)BN_ALLOC(sizeof(uint32_t) * m_capacity);
		memset(digits, 0, m_size << 2);

		uint32_t ls = 0;
		for (uint32_t* d = b.digits + b.m_size; d != b.digits;) a.addmul(*(--d), ls++, *this);

		pack();
	}
}

inline void bignatural::addmul(uint32_t b, uint32_t left_shifts, bignatural& result) const {
	uint64_t carry = 0;
	uint32_t* dp = digits + m_size;
	uint32_t* rp = result.digits + result.m_size - left_shifts;
	do {
		uint64_t pm = ((uint64_t)(*(--dp))) * b + carry + (*(--rp));
		*rp = (uint32_t)pm;
		carry = pm >> 32;
	} while (dp != digits);

	*(--rp) = (uint32_t)carry;
}

inline void bignatural::push_back(uint32_t b) {
	if (m_size == m_capacity) increaseCapacity((m_capacity | 1) << 2);
	digits[m_size++] = b;
}

inline void bignatural::push_bit_back(uint32_t b) {
	if (m_size) {
		leftShift(1);
		back() |= b;
	}
	else if (b) push_back(1);
}

inline void bignatural::pop_front() {
	uint32_t* d = digits, * de = digits + m_size;
	while (++d != de) *(d - 1) = *d;
	pop_back();
}

inline void bignatural::pack() {
	uint32_t i = 0;
	while (i < m_size && digits[i] == 0) i++;

	if (i) {
		uint32_t* dold = digits + i;
		uint32_t* dnew = digits;
		uint32_t* dend = digits + m_size;
		while (dold < dend) *dnew++ = *dold++;
		m_size -= i;
	}
}


inline bignatural& bignatural::operator+=(const uint32_t b) {
	if (b != 0) {
		if (m_size) {
			uint64_t pm = ((uint64_t)digits[(int)(m_size - 1)]) + b;
			digits[(int)m_size - 1] = (uint32_t)(pm);
			uint64_t carry = pm >> 32;
			for (uint32_t i = m_size - 1; carry && i > 0; i--) {
				pm = ((uint64_t)digits[(int)(i - 1)]) + carry;
				digits[(int)i - 1] = (uint32_t)(pm);
				carry = pm >> 32;
			}
			if (carry) addOneMostSignificantDigit((uint32_t)carry);
		}
		else init(b);
	}

	return *this;
}

inline bignatural& bignatural::operator+=(const uint64_t b) {
	if (b <= UINT32_MAX) return operator+=((uint32_t)b);
	else return operator+=(bignatural(b));
}

inline bignatural& bignatural::operator*=(const uint32_t b) {
	if (m_size) {
		if (b == 0) return operator=(0U);
		uint32_t* da = digits + m_size;
		uint64_t pm = ((uint64_t)(*(--da))) * b;
		*da = (uint32_t)(pm);
		uint64_t carry = pm >> 32;
		while (da != digits) {
			pm = ((uint64_t)(*(--da))) * b + carry;
			*da = (uint32_t)(pm);
			carry = pm >> 32;
		}

		if (carry) addOneMostSignificantDigit((uint32_t)carry);
	}

	return *this;
}

inline bignatural& bignatural::operator*=(const uint64_t b) {
	if (b < UINT32_MAX) return operator*=((uint32_t)b);
	else return operator*=(bignatural(b));
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// FUNCTION IMPLEMENTATIONS - BIGFLOAT
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

inline bigfloat::bigfloat() : exponent(0), sign(0) {}

inline bigfloat::bigfloat(const bignatural& m, int32_t e, int32_t s) : mantissa(m), exponent(e), sign(s) {}

inline int32_t bigfloat::log2() const { return exponent + ((int32_t)mantissa.getNumSignificantBits()) - 1; }

inline void bigfloat::increaseMantissa() { mantissa += 1U; pack(); }

inline bool bigfloat::operator!=(const bigfloat& b) const { return (operator-(b).sign != 0); }

inline void bigfloat::invert() { sign = -sign; }

inline int bigfloat::sgn() const { return sign; }

inline const bignatural& bigfloat::getMantissa() const { return mantissa; }

inline int32_t bigfloat::getExponent() const { return exponent; }

inline void bigfloat::leftShift(uint32_t n) {
	mantissa <<= n;
	exponent -= n;
}

inline void bigfloat::rightShift(uint32_t n) {
	mantissa >>= n;
	exponent += n;
}

inline bigfloat::bigfloat(const double d) {
	sign = (d > 0) - (d < 0);

	if (sign) {
		uint64_t dn = *((uint64_t*)(&d));
		const uint64_t m = (dn & 0x000fffffffffffff) + 0x0010000000000000;
		mantissa.push_back(m >> 32);
		mantissa.push_back((uint32_t)m);
		dn <<= 1;
		dn >>= 53;
		exponent = ((int32_t)dn) - 1075; // Exp

		pack();
	}
	else exponent = 0;
}

inline double bigfloat::get_d() const {
	uint64_t dn = 0;
	if (mantissa.empty()) return 0.0;

	uint64_t m;
	int32_t e;
	uint32_t shft;

	if (mantissa.size() == 1) {
		m = ((uint64_t)mantissa[0]);
		shft = mantissa.countLeadingZeroes() + 21;
		m <<= shft;
		e = exponent - (int32_t)shft;
	}
	else {
		m = (((uint64_t)mantissa[0]) << 32) | ((uint64_t)mantissa[1]);
		e = exponent + 32 * ((int32_t)mantissa.size() - 2);
		shft = mantissa.countLeadingZeroes();

		if (shft < 11) {
			m >>= (11 - shft);
			e += (11 - shft);
		}
		if (shft > 11) {
			m <<= (shft - 11);
			e -= (shft - 11);
			if (mantissa.size() > 2) m |= (mantissa[2] >> (43 - shft));
		}
	}
	m &= (~0x0010000000000000); // Remove implicit digit
	e += 52;

	if (e < (-1022)) return 0.0;
	if (e > 1023) return ((double)sign) * INFINITY;

	if (sign < 0) dn |= 0x8000000000000000; // Set sign
	dn |= (((uint64_t)(e + 1023)) << 52); // Set exponent
	dn |= m; // Set mantissa

	return *((double*)(&dn));
}

inline bigfloat bigfloat::operator+(const bigfloat& b) const {
	if (mantissa.empty()) return b;
	if (b.mantissa.empty()) return *this;

	if (exponent == b.exponent) {
		bigfloat result;

		if (sign == b.sign) {
			result.mantissa.toSum(mantissa, b.mantissa);
			result.sign = sign;
		}
		else if (b.mantissa >= mantissa) {
			result.mantissa.toDiff(b.mantissa, mantissa);
			result.sign = b.sign;
		}
		else {
			result.mantissa.toDiff(mantissa, b.mantissa);
			result.sign = sign;
		}

		result.exponent = exponent;
		result.pack();
		return result;
	}
	else if (exponent > b.exponent) {
		bigfloat op(*this);
		op.leftShift((uint32_t)(exponent - b.exponent));
		return op + b;
	}
	else { // exponent < b.exponent
		bigfloat op(b);
		op.leftShift((uint32_t)(b.exponent - exponent));
		return op + *this;
	}
}

inline bigfloat bigfloat::operator-(const bigfloat& b) const {
	if (mantissa.empty()) return b.inverse();
	if (b.mantissa.empty()) return *this;

	if (exponent == b.exponent) {
		bigfloat result;

		if (sign != b.sign) {
			result.mantissa.toSum(mantissa, b.mantissa);
			result.sign = sign;
		}
		else if (b.mantissa >= mantissa) {
			result.mantissa.toDiff(b.mantissa, mantissa);
			result.sign = -sign;
		}
		else {
			result.mantissa.toDiff(mantissa, b.mantissa);
			result.sign = sign;
		}

		result.exponent = exponent;
		result.pack();
		return result;
	}
	else if (exponent > b.exponent) {
		bigfloat op(*this);
		op.leftShift((uint32_t)(exponent - b.exponent));
		return op - b;
	}
	else { // exponent < b.exponent
		bigfloat op(b);
		op.leftShift((uint32_t)(b.exponent - exponent));
		return *this - op;
	}
}


inline bigfloat bigfloat::operator*(const bigfloat& b) const {
	if (mantissa.empty() || b.mantissa.empty()) return 0;

	bigfloat result;
	result.mantissa.toProd(mantissa, b.mantissa);
	result.sign = sign * b.sign;
	result.exponent = exponent + b.exponent;
	result.pack();
	return result;
}

inline std::string bigfloat::get_str() const {
	std::string s;
	if (sign == 0) s += "0";
	if (sign < 0) s += "-";
	s += mantissa.get_str();
	s += " * 2^";
	s += std::to_string(exponent);
	return s;
}

inline void bigfloat::pack() {
	if (mantissa.empty()) {
		sign = exponent = 0;
		return;
	}

	while (mantissa.back() == 0) {
		mantissa.pop_back();
		exponent += 32;
	}

	const uint32_t s = mantissa.countEndingZeroesLSL();
	if (s) {
		auto dp = mantissa.digits, de = mantissa.digits + mantissa.m_size - 1;
		const uint32_t ts = 32 - s;
		while (dp != de) {
			*de >>= s;
			*de |= ((*(de - 1)) << ts);
			de--;
		}
		*de >>= s;
		exponent += s;
	}

	mantissa.pack();
}

inline bigfloat bigfloat::sqrt(uint32_t prec_bits) const {
	assert(sign >= 0);
	bigfloat a = *this;

	// Need twice precision here, because the sqrt will halve it
	uint32_t shf = prec_bits * 2;

	// Possibly shift the mantissa to get the precision required
	uint32_t msb = a.mantissa.getNumSignificantBits();
	if (shf > msb) a.leftShift(shf - msb);
	else a.rightShift(msb - shf);

	// If exponent is not divisible by 2, add 1 and shift
	if (a.exponent & 1) a.leftShift(1);

	// Now compute the sqrt
	// exponent is halved and mantissa is sqrt'ed
	a.exponent >>= 1;
	a.mantissa = a.mantissa.sqrt();
	a.pack();

	return a;
}

inline bigfloat bigfloat::inverse() const {
	bigfloat r = *this;
	r.invert();
	return r;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// FUNCTION IMPLEMENTATIONS - BIGRATIONAL
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

inline void bigrational::invert() {
	assert(sign != 0);
	std::swap(numerator, denominator);
}

inline bigrational bigrational::inverse() const { bigrational r = *this; r.invert(); return r; }

inline void bigrational::negate() { sign = -sign; }

inline bigrational bigrational::negation() const { bigrational r = *this; r.negate(); return r; }

inline bigrational bigrational::operator-(const bigrational& r) const { return operator+(r.negation()); }

inline bigrational bigrational::operator*(const bigrational& r) const {
	if (sign == 0 || r.sign == 0) return bigrational();
	else return bigrational(numerator * r.numerator, denominator * r.denominator, sign * r.sign);
}

inline bigrational bigrational::operator/(const bigrational& r) const {
	assert(r.sign != 0);
	return operator*(r.inverse());
}

inline bool bigrational::operator==(const bigrational& r) const {
	return (sign == r.sign && numerator == r.numerator && denominator == r.denominator);
}

inline bool bigrational::operator!=(const bigrational& r) const {
	return (sign != r.sign || numerator != r.numerator || denominator != r.denominator);
}

inline bool bigrational::operator>(const bigrational& r) const {
	return (sign > r.sign || (sign > 0 && r.sign > 0 && hasGreaterModule(r)) || (sign < 0 && r.sign < 0 && r.hasGreaterModule(*this)));
}

inline bool bigrational::operator>=(const bigrational& r) const {
	return (sign > r.sign || (sign > 0 && r.sign > 0 && hasGrtrOrEqModule(r)) || (sign < 0 && r.sign < 0 && r.hasGrtrOrEqModule(*this)));
}

inline bool bigrational::operator<(const bigrational& r) const {
	return (sign < r.sign || (sign < 0 && r.sign < 0 && hasGreaterModule(r)) || (sign > 0 && r.sign > 0 && r.hasGreaterModule(*this)));
}

inline bool bigrational::operator<=(const bigrational& r) const {
	return (sign < r.sign || (sign < 0 && r.sign < 0 && hasGrtrOrEqModule(r)) || (sign > 0 && r.sign > 0 && r.hasGrtrOrEqModule(*this)));
}

inline bool bigrational::hasGreaterModule(const bigrational& r) const {
	return numerator * r.denominator > r.numerator * denominator;
}

inline bool bigrational::hasGrtrOrEqModule(const bigrational& r) const {
	return numerator * r.denominator >= r.numerator * denominator;
}

inline const bignatural& bigrational::get_num() const { return numerator; }
inline const bignatural& bigrational::get_den() const { return denominator; }

inline int32_t bigrational::sgn() const { return sign; }

inline bigrational::bigrational(const bigfloat& f) {
	if (f.sgn() == 0) sign = 0;
	else {
		sign = f.sgn();
		numerator = f.getMantissa();
		denominator = 1;
		int32_t e = f.getExponent();
		if (e >= 0) numerator.leftShift((uint32_t)e);
		else denominator.leftShift((uint32_t)(-e));
	}
}

inline void bigrational::init(FILE* fp) {
	sign = 1;
	int c;
	while (isspace(c = fgetc(fp)));
	if (c == '-') sign = -1;
	else if (c != EOF) ungetc(c, fp);

	numerator = bignatural(fp);

	c = fgetc(fp);
	if (c == '/') {
		denominator = bignatural(fp);
		if (numerator.empty()) {
			sign = 0;
			denominator = 0;
		}
	}
	else if (c == '.') {
		uint64_t num_digits = (uint64_t)ftell(fp);
		bignatural decimal_part = bignatural(fp);
		num_digits = ftell(fp) - num_digits; // Number of decimal digits read
		denominator = 1U;
		while (num_digits--) denominator *= 10U;
		numerator *= denominator;
		numerator += decimal_part;
	}
	else {
		if (c != EOF) ungetc(c, fp);
		denominator = 1;
	}
}

inline void bigrational::compress() {
	const uint32_t nez = numerator.countEndingZeroes();
	const uint32_t dez = denominator.countEndingZeroes();
	const uint32_t s = std::min(nez, dez);
	numerator >>= s;
	denominator >>= s;
}

inline void bigrational::canonicalize() {
	assert(denominator.size());
	if (sign) {
		if (numerator.empty()) {
			numerator = denominator = 0;
			sign = 0;
		}
		else {
			compress();
			bignatural r;
			const bignatural gcd = numerator.GCD(denominator);
			numerator = numerator.divide_by(gcd, r);
			denominator = denominator.divide_by(gcd, r);
		}
	}
}

inline bigrational bigrational::operator+(const bigrational& r) const {
	if (sign == 0) return r;
	else if (r.sign == 0) return *this;
	else {
		//bignatural rm;
		//const bignatural gcd = denominator.GCD(r.denominator);
		//const bignatural den3 = (denominator * r.denominator).divide_by(gcd, rm);
		//const bignatural left_den = den3.divide_by(denominator, rm);
		//const bignatural right_den = den3.divide_by(r.denominator, rm);
		//const bignatural left_num = numerator * left_den;
		//const bignatural right_num = r.numerator * right_den;
		//if (sign > 0 && r.sign > 0)	return bigrational(left_num + right_num, den3, 1);
		//else if (sign < 0 && r.sign < 0) return bigrational(left_num + right_num, den3, -1);
		//else if (sign > 0 && r.sign < 0) {
		//	if (left_num >= right_num) return bigrational(left_num - right_num, den3, 1);
		//	else return bigrational(right_num - left_num, den3, -1);
		//}
		//else { // if (sign < 0 && r.sign > 0)
		//	if (left_num >= right_num) return bigrational(left_num - right_num, den3, -1);
		//	else return bigrational(right_num - left_num, den3, 1);
		//}

		const bignatural left_num = numerator * r.denominator;
		const bignatural right_num = r.numerator * denominator;
		if (sign > 0 && r.sign > 0)	return bigrational(left_num + right_num, denominator * r.denominator, 1);
		else if (sign < 0 && r.sign < 0) return bigrational(left_num + right_num, denominator * r.denominator, -1);
		else if (sign > 0 && r.sign < 0) {
			if (left_num >= right_num) return bigrational(left_num - right_num, denominator * r.denominator, 1);
			else return bigrational(right_num - left_num, denominator * r.denominator, -1);
		}
		else { // if (sign < 0 && r.sign > 0)
			if (left_num >= right_num) return bigrational(left_num - right_num, denominator * r.denominator, -1);
			else return bigrational(right_num - left_num, denominator * r.denominator, 1);
		}
	}
}

inline double bigrational::get_d() const {
	if (sign == 0) return 0.0;

	bignatural num = numerator;
	bignatural den = denominator;
	int32_t E = (int32_t)num.getNumSignificantBits() - (int32_t)den.getNumSignificantBits();
	if (E > 0) { den.leftShift((uint32_t)E); if (den > num) { E--; den >>= 1; } }
	else if (E <= 0) { num.leftShift((uint32_t)(-E)); if (den > num) { E--; num.leftShift(1); } }

	if (E > 1023) return INFINITY;
	else if (E < -1022) return -INFINITY;

	uint64_t signbit = sign < 0 ? ((uint64_t)1) << 63 : 0;
	uint64_t exponent = ((uint64_t)(1023 + E)) << 52;
	uint64_t mantissa = 0;

	for (int i = 0; i < 53; i++) {
		mantissa <<= 1;
		if (num >= den) {
			mantissa |= 1;
			num = num - den;
		}
		num <<= 1;
	}
	mantissa &= (~(((uint64_t)1) << 52));
	mantissa |= exponent;
	mantissa |= signbit;

	void* ptr = &mantissa;

	return *((double*)ptr);
}

inline bigfloat bigrational::get_bigfloat(uint32_t num_significant_bits) const {
	if (sign == 0) return 0.0;

	bignatural num = numerator;
	bignatural den = denominator;
	int32_t E = (int32_t)num.getNumSignificantBits() - (int32_t)den.getNumSignificantBits();
	if (E > 0) { den.leftShift((uint32_t)E); if (den > num) { E--; den >>= 1; } }
	else if (E <= 0) { num.leftShift((uint32_t)(-E)); if (den > num) { E--; num.leftShift(1); } }

	bignatural mantissa;

	for (uint32_t i = 0; i <= num_significant_bits; i++) {
		mantissa <<= 1;
		if (num >= den) {
			mantissa |= 1;
			num = num - den;
		}
		num <<= 1;
	}

	return bigfloat(mantissa, E - (int32_t)num_significant_bits, sign);
}

inline std::string bigrational::get_dec_str() const {
	std::string st;
	if (sign < 0) st = "-";
	else if (sign == 0) return "0";
	st += numerator.get_dec_str();
	const std::string dens = denominator.get_dec_str();
	if (dens != "1") {
		st += "/";
		st += dens;
	}
	return st;
}

inline std::string bigrational::get_str() const {
	std::string st;
	if (sign < 0) st = "-";
	else if (sign == 0) return "0";
	st += numerator.get_str();
	st += "/";
	st += denominator.get_str();
	return st;
}

#endif //USE_GNU_GMP_CLASSES
