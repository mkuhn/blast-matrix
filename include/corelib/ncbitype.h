#ifndef CORELIB___NCBITYPE__H
#define CORELIB___NCBITYPE__H

/*  $Id: ncbitype.h 164289 2009-06-24 21:18:16Z vakatov $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Denis Vakatov
 *
 *
 */

/**
 * @file ncbitype.h
 *
 * Defines NCBI C/C++ fixed-size types.
 *
 *   -  Char, Uchar
 *   -  Int1, Uint1
 *   -  Int2, Uint2
 *   -  Int4, Uint4
 *   -  Int8, Uint8
 *   -  Ncbi_BigScalar
 *   -  Macros for constant values definition.
 *
 */

/** @addtogroup Portability
 *
 * @{
 */

#include <ncbiconf.h>

#ifdef HAVE_INTTYPES_H
#  include <inttypes.h>
#endif


/* For the sake of backward compatibility only */
#include <common/ncbiconf_impl.h>


/* Char, Uchar, Int[1,2,4], Uint[1,2,4]
 */

#if (SIZEOF_CHAR != 1)
#  error "Unsupported size of char(must be 1 byte)"
#endif
#if (SIZEOF_SHORT != 2)
#  error "Unsupported size of short int(must be 2 bytes)"
#endif
#if (SIZEOF_INT != 4)
#  error "Unsupported size of int(must be 4 bytes)"
#endif


typedef          char  Char;    /**< Alias for char */
typedef signed   char  Schar;   /**< Alias for signed char */
typedef unsigned char  Uchar;   /**< Alias for unsigned char */
typedef signed   char  Int1;    /**< Alias for signed char */
typedef unsigned char  Uint1;   /**< Alias for unsigned char */
typedef signed   short Int2;    /**< Alias for signed short */
typedef unsigned short Uint2;   /**< Alias for unsigned short */
typedef signed   int   Int4;    /**< Alias for signed int */
typedef unsigned int   Uint4;   /**< Alias for unsigned int */


/* Int8, Uint8
 */

#if   (SIZEOF_LONG == 8)
#  define NCBI_INT8_TYPE         long
#  define NCBI_INT8_IS_LONG      1
#elif (SIZEOF_LONG_LONG == 8)
#  define NCBI_INT8_TYPE         long long
#  define NCBI_INT8_IS_LONG_LONG 1
#elif (SIZEOF___INT64 == 8)
#  define NCBI_INT8_TYPE         __int64
#  define NCBI_INT8_IS_INT64     1
/** @deprecated  Use NCBI_INT8_IS_INT64 instead */
#  define NCBI_USE_INT64         1
#else
#  error "This platform does not support 8-byte integer"
#endif

/** Signed 8 byte sized integer */
typedef signed   NCBI_INT8_TYPE Int8;    

/** Unsigned 8 byte sized integer */
typedef unsigned NCBI_INT8_TYPE Uint8;


/* BigScalar
 */

#define NCBI_BIG_TYPE NCBI_INT8_TYPE
#define SIZEOF_NCBI_BIG 8
#if (SIZEOF_LONG_DOUBLE > SIZEOF_NCBI_BIG)
#  undef  NCBI_BIG_TYPE
#  undef  SIZEOF_NCBI_BIG
#  define NCBI_BIG_TYPE   long double 
#  define SIZEOF_NCBI_BIG SIZEOF_LONG_DOUBLE
#endif
#if (SIZEOF_DOUBLE > SIZEOF_NCBI_BIG)
#  undef  NCBI_BIG_TYPE
#  undef  SIZEOF_NCBI_BIG
#  define NCBI_BIG_TYPE   double
#  define SIZEOF_NCBI_BIG SIZEOF_DOUBLE
#endif
#if (SIZEOF_VOIDP > SIZEOF_NCBI_BIG)
#  undef  NCBI_BIG_TYPE
#  undef  SIZEOF_NCBI_BIG
#  define NCBI_BIG_TYPE   void*
#  define SIZEOF_NCBI_BIG SIZEOF_VOIDP
#endif

/**
 * Define large scalar type.
 *
 * This is platform dependent. It could be an Int8, long double, double
 * or void*.
 */
typedef NCBI_BIG_TYPE Ncbi_BigScalar;


#ifndef HAVE_INTPTR_T
#  if SIZEOF_INT == SIZEOF_VOIDP
typedef int intptr_t;
#  elif SIZEOF_LONG == SIZEOF_VOIDP
typedef long intptr_t;
#  elif SIZEOF_LONG_LONG == SIZEOF_VOIDP
typedef long long intptr_t;
#  else
#    error No integer type is the same size as a pointer!
#  endif
#endif

#ifndef HAVE_UINTPTR_T
#  if SIZEOF_INT == SIZEOF_VOIDP
typedef unsigned int uintptr_t;
#  elif SIZEOF_LONG == SIZEOF_VOIDP
typedef unsigned long uintptr_t;
#  elif SIZEOF_LONG_LONG == SIZEOF_VOIDP
typedef unsigned long long uintptr_t;
#  else
#    error No integer type is the same size as a pointer!
#  endif
#endif


/* Macros for constant values definition 
 */

#if (SIZEOF_LONG == 8)
#  define NCBI_CONST_INT8(v)   v##L
#  define NCBI_CONST_UINT8(v)  v##UL
#elif (SIZEOF_LONG_LONG == 8)
#  define NCBI_CONST_INT8(v)   v##LL
#  define NCBI_CONST_UINT8(v)  v##ULL
#elif defined(NCBI_USE_INT64)
#  define NCBI_CONST_INT8(v)   v##i64
#  define NCBI_CONST_UINT8(v)  v##ui64
#else
#  define NCBI_CONST_INT8(v)   v
#  define NCBI_CONST_UINT8(v)  v
#endif


/* Undef auxiliaries
 */

#undef SIZEOF_NCBI_BIG
#undef NCBI_BIG_TYPE
#undef NCBI_INT8_TYPE


#endif  /* CORELIB___NCBITYPE__H */


/* @} */
