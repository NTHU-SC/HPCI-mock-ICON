/**
 * @file xt_arithmetic_util.h
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef XT_ARITHMETIC_UTIL_H
#define XT_ARITHMETIC_UTIL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <limits.h>
#include <stddef.h>

#if HAVE_DECL___LZCNT || HAVE_DECL___LZCNT64 || HAVE_DECL___LZCNT16
#  include <intrin.h>
#endif

#include "xt/xt_core.h"

enum {
  xt_int_bits = sizeof (Xt_int) * CHAR_BIT,
};

/* simple operations on Xt_int */
/**
 * @return -1 if x < 0, 1 otherwise
 */
static inline Xt_int
Xt_isign(Xt_int x)
{
#if (-1 >> 1) == -1
  return ((x >> (sizeof (Xt_int) * CHAR_BIT - 1)) |
          (Xt_int)((Xt_uint)(~x) >> (sizeof (Xt_int) * CHAR_BIT - 1)));
#else
  return (Xt_int)((x >= 0) - (x < 0));
#endif
}

/**
 * @return -1 if x < 0, 1 otherwise
 */
static inline int
isign(int x)
{
#if (-1 >> 1) == -1
  return ((x >> (sizeof (x) * CHAR_BIT - 1)) |
          (int)((unsigned int)(~x) >> (sizeof (x) * CHAR_BIT - 1)));
#else
  return (x >= 0) - (x < 0);
#endif
}

/**
 * @return ~0 if x < 0, 0 otherwise
 */
static inline int
isign_mask(int x)
{
#if (-1 >> 1) == -1
  return x >> (sizeof (int) * CHAR_BIT - 1);
#else
#warning Unusual behaviour of shift operator detected.
  return (x < 0) * ~0;
#endif
}

/**
 * @return ~0 if x < 0, 0 otherwise
 */
static inline MPI_Aint
asign_mask(MPI_Aint x)
{
#if (-1 >> 1) == -1
  return x >> (sizeof (MPI_Aint) * CHAR_BIT - 1);
#else
#warning Unusual behaviour of shift operator detected.
  return (x < 0) * ~(MPI_Aint)0;
#endif
}


/**
 * @return ~(Xt_int)0 if x < 0, 0 otherwise
 */
static inline Xt_int
Xt_isign_mask(Xt_int x)
{
#if (-1 >> 1) == -1
  return x >> (sizeof (x) * CHAR_BIT - 1);
#else
#warning Unusual behaviour of shift operator detected.
  return (x < 0) * ~(Xt_int)0;
#endif
}

/**
 * @return -1 if x < 0, 1 otherwise
 */
static inline long long
llsign(long long x)
{
#if (-1 >> 1) == -1
  return ((x >> (sizeof (x) * CHAR_BIT - 1)) |
          (long long)((unsigned long long)(~x) >> (sizeof (x) * CHAR_BIT - 1)));
#else
  return (x >= 0) - (x < 0);
#endif
}

/**
 * @return ~(long long)0 if x < 0, 0 otherwise
 */
static inline long long
llsign_mask(long long x)
{
#if (-1 >> 1) == -1
  return x >> (sizeof (x) * CHAR_BIT - 1);
#else
#warning Unusual behaviour of shift operator detected.
  return (x < 0) * ~0LL;
#endif
}

/**
 * @return MIN(a, b)
 */
static inline int
imin(int a, int b)
{
  return a <= b ? a : b;
}

/* return number of leading zeroes in an Xt_uint */
#ifdef XT_INT_CLZ
#define xinlz(v) XT_INT_CLZ(v)
#else
static inline int
xinlz(Xt_uint v)
{
  int c = 0;
#if SIZEOF_XT_INT * CHAR_BIT == 64
  if (v <= UINT64_C(0x00000000ffffffff)) {c += 32; v <<= 32;}
  if (v <= UINT64_C(0x0000ffffffffffff)) {c += 16; v <<= 16;}
  if (v <= UINT64_C(0x00ffffffffffffff)) {c +=  8; v <<=  8;}
  if (v <= UINT64_C(0x0fffffffffffffff)) {c +=  4; v <<=  4;}
  if (v <= UINT64_C(0x3fffffffffffffff)) {c +=  2; v <<=  2;}
  if (v <= UINT64_C(0x7fffffffffffffff)) {c +=  1;}
#elif SIZEOF_XT_INT * CHAR_BIT == 32
  if (v <= 0x0000ffffUL) {c += 16; v <<= 16;}
  if (v <= 0x00ffffffUL) {c +=  8; v <<=  8;}
  if (v <= 0x0fffffffUL) {c +=  4; v <<=  4;}
  if (v <= 0x3fffffffUL) {c +=  2; v <<=  2;}
  if (v <= 0x7fffffffUL) {c +=  1;}
#elif SIZEOF_XT_INT * CHAR_BIT == 16
  if (v <= 0x00ffU) {c +=  8; v <<=  8;}
  if (v <= 0x0fffU) {c +=  4; v <<=  4;}
  if (v <= 0x3fffU) {c +=  2; v <<=  2;}
  if (v <= 0x7fffU) {c +=  1;}
#else
#error "Unexpected size of Xt_int.\n"
#endif
  return c;
}
#endif

/* return number of trailing zeroes in an Xt_uint */
#ifdef XT_INT_CTZ
#define xintz(v) XT_INT_CTZ(v)
#else
static inline int
xintz(Xt_uint v)
{
  Xt_uint lc = (Xt_uint)((v ^ (v - 1)) >> 1);
  return lc ? xt_int_bits - xinlz(lc) : 0;
}
#endif

#ifdef HAVE_ASM_BSR
static inline size_t
next_2_pow(size_t v)
{
  enum {
    size_t_bits = sizeof (size_t) * CHAR_BIT,
  };
  size_t r;
  if (v <= 1) {
    r = 1;
  } else {
    size_t ms1bpos;
#if SIZEOF_LONG == 8
    __asm__ ("bsrq %1, %0" : "=r" (ms1bpos) : "r" (v-1));
#elif SIZEOF_LONG == 4
    __asm__ ("bsrl %1, %0" : "=r" (ms1bpos) : "r" (v-1));
#else
#error "Unexpected size of size_t!"
#endif
    r = (size_t)1 << (ms1bpos+1);
  }
  return r;
}
#elif HAVE_DECL___BUILTIN_CLZL
#define clzl(v) (__builtin_clzl(v))
#elif HAVE_DECL___LZCNT && SIZEOF_LONG == SIZEOF_INT
#define clzl(v) ((int)(__lzcnt(v)))
#elif HAVE_DECL___LZCNT64 && SIZEOF_LONG == 8 && CHAR_BIT == 8
#define clzl(v) ((int)(__lzcnt64(v)))
#else
static inline size_t
next_2_pow(size_t v)
{
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
#if SIZEOF_SIZE_T * CHAR_BIT == 64
  v |= v >> 32;
#endif
  return v+1;
}
#endif
#ifdef clzl
static inline size_t
next_2_pow(size_t v)
{
  enum {
    size_t_bits = sizeof (size_t) * CHAR_BIT,
  };
  size_t r;
  if (v == 0) {
    r = 1;
  } else {
    r = clzl(v-1);
    r = (size_t)1 << (size_t_bits - r);
  }
  return r;
}

#endif

static inline Xt_int
Xt_doz(Xt_int a, Xt_int b)
{
  return (Xt_int)((a - b) & -(a >= b));
}


// For signed integers, a similar method follows.
//
// Given c > 1 and odd, compute m such that (c * m) mod 2^n == 1
// Then if c divides x (x%c ==0), the quotient is given by q = x/c == x*m mod 2^n
//
// x can range from ⎡-2^(n-1)/c⎤ * c, ... -c, 0, c, ...  ⎣(2^(n-1) - 1)/c⎦ * c
// Thus, x*m mod 2^n is ⎡-2^(n-1)/c⎤, ... -2, -1, 0, 1, 2, ... ⎣(2^(n-1) - 1)/c⎦
//
// So, x is a multiple of c if and only if:
// ⎡-2^(n-1)/c⎤ <= x*m mod 2^n <= ⎣(2^(n-1) - 1)/c⎦
//
// Since c > 1 and odd, this can be simplified by
// ⎡-2^(n-1)/c⎤ == ⎡(-2^(n-1) + 1)/c⎤ == -⎣(2^(n-1) - 1)/c⎦
//
// -⎣(2^(n-1) - 1)/c⎦ <= x*m mod 2^n <= ⎣(2^(n-1) - 1)/c⎦
//
// To extend this to even integers, consider c = d0 * 2^k where d0 is odd.
// We can test whether x is divisible by both d0 and 2^k.
//
// Let m be such that (d0 * m) mod 2^n == 1.
// Let q = x*m mod 2^n. Then c divides x if:
//
// -⎣(2^(n-1) - 1)/d0⎦ <= q <= ⎣(2^(n-1) - 1)/d0⎦ and q ends in at least k 0-bits
//
// To transform this to a single comparison, we use the following theorem (ZRS in Hacker's Delight).
//
// For a >= 0 the following conditions are equivalent:
// 1) -a <= x <= a and x ends in at least k 0-bits
// 2) RotRight(x+a', k) <= ⎣2a'/2^k⎦
//
// Where a' = a & -2^k (a with its right k bits set to zero)
//
// To see that 1 & 2 are equivalent, note that -a <= x <= a is equivalent to
// -a' <= x <= a' if and only if x ends in at least k 0-bits.  Adding -a' to each side gives,
// 0 <= x + a' <= 2a' and x + a' ends in at least k 0-bits if and only if x does since a' has
// k 0-bits by definition.  We can use theorem ZRU above with x -> x + a' and a -> 2a' giving 1) == 2).
//
// Let m be such that (d0 * m) mod 2^n == 1.
// Let q = x*m mod 2^n.
// Let a' = ⎣(2^(n-1) - 1)/d0⎦ & -2^k
//
// Then the divisibility test is:
//
// RotRight(q+a', k) <= ⎣2a'/2^k⎦
//
// Note that the calculation is performed using unsigned integers.
// Since a' can have n-1 bits, 2a' may have n bits and there is no
// risk of overflow.

#if defined XT_USE_FAST_DIVISIBLE_TEST \
  || (defined __ICC && defined __OPTIMIZE__)    \
  || (defined __GNUC__ && defined __OPTIMIZE__) \
  || (defined __PGI)
#undef XT_USE_FAST_DIVISIBLE_TEST
#define XT_USE_FAST_DIVISIBLE_TEST

struct xt_fast_div_coeff {
  Xt_uint bs;  // trailingZeros(c)
  Xt_uint minv; // minv * (c>>bs) mod 2^n == 1 multiplicative inverse of
                // odd portion modulo 2^n
  Xt_uint a;    // ⎣(2^(n-1) - 1)/ (c>>bs)⎦ & -(1<<bs) additive constant
  Xt_uint max;  // ⎣(2 a) / (1<<bs)⎦ max value to for divisibility
};

static inline struct xt_fast_div_coeff
get_fast_div_coeff(Xt_uint d)
{
  assert(d > 0);
  // the initial values work for d == 1
  struct xt_fast_div_coeff coeff = { .bs = 0, .minv = 0, .a = 0, .max = 0 };
  if (d > 1) {
    int bs = xintz(d);
    Xt_uint d0 = (Xt_uint)(d >> bs); // the odd portion of the divisor

    Xt_uint mask = (Xt_uint)(~(Xt_uint)0);

    // Calculate the multiplicative inverse via Newton's method.
    // Quadratic convergence doubles the number of correct bits per iteration.
    Xt_uint m = d0;            // initial guess correct to 3-bits d0*d0 mod 8 == 1
    m = (Xt_uint)(m * (2 - m*d0)); // 6-bits
    m = (Xt_uint)(m * (2 - m*d0)); // 12-bits
    m = (Xt_uint)(m * (2 - m*d0)); // 24-bits
#if SIZEOF_XT_INT * CHAR_BIT > 24
    m = (Xt_uint)(m * (2 - m*d0)); // 48-bits
#if SIZEOF_XT_INT * CHAR_BIT > 48
    m = (Xt_uint)(m * (2 - m*d0)); // 96-bits >= 64-bits
#endif
#endif
    Xt_uint a = (Xt_uint)(((mask >> 1) / d0) & (Xt_uint)(-(1 << bs)));
    Xt_uint max = (Xt_uint)((2 * a) >> bs);

    coeff.bs = (Xt_uint)bs;
    coeff.minv = m;
    coeff.a = a;
    coeff.max = max;
  }
  return coeff;
}

static inline int
fast_divisible(struct xt_fast_div_coeff coeff, Xt_int i)
{
  Xt_uint k = coeff.bs;
  Xt_uint mul = (Xt_uint)((Xt_uint)i * coeff.minv + coeff.a);
  Xt_uint rot = (Xt_uint)((mul >> k)
                          | (mul<<((xt_int_bits-k)&(xt_int_bits-1))));
  return rot <= coeff.max;
}

#define is_divisible(divisor, coeff, dividend) fast_divisible(coeff, dividend)
#else
#define is_divisible(divisor, coeff, dividend) (divisor==0 || (dividend)%(divisor)==0)
#endif

#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
