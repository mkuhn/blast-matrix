#ifndef CONNECT___NCBI_ANSI_EXT__H
#define CONNECT___NCBI_ANSI_EXT__H

/* $Id: ncbi_ansi_ext.h 192302 2010-05-22 14:54:54Z lavr $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   Non-ANSI, yet widely used functions
 *
 */

#include <connect/connect_export.h>
#include "ncbi_config.h"
#include <stddef.h>
#include <string.h>


#ifdef __cplusplus
extern "C" {
#endif


#ifndef HAVE_STRDUP

#  ifdef   strdup
#    undef strdup
#  endif
#  define  strdup      NCBI_strdup

/* Create a copy of string "str".
 * Return an identical malloc'ed string, which must be explicitly freed 
 * by free() when no longer needed.
 */
extern NCBI_XCONNECT_EXPORT
char* strdup(const char* str);

#endif /*HAVE_STRDUP*/


#ifndef HAVE_STRNDUP

#  ifdef   strndup
#    undef strndup
#  endif
#  define  strndup     NCBI_strndup

/* Create a copy of up to "n" first characters of string "str".
 * Return a malloc'ed and '\0'-terminated string, which must be
 * explicitly freed by free() when no longer needed.
 */
extern NCBI_XCONNECT_EXPORT
char* strndup(const char* str, size_t n);

#endif /*HAVE_STRNDUP*/


#ifndef HAVE_STRCASECMP

#  ifdef   strcasecmp
#    undef strcasecmp
#    undef strncasecmp
#  endif
#  define  strcasecmp  NCBI_strcasecmp
#  define  strncasecmp NCBI_strncasecmp

/* Compare "s1" and "s2", ignoring case.
 * Return less than, equal to or greater than zero if
 * "s1" is lexicographically less than, equal to or greater than "s2".
 */
extern NCBI_XCONNECT_EXPORT
int strcasecmp(const char* s1, const char* s2);

/* Compare not more than "n" characters of "s1" and "s2", ignoring case.
 * Return less than, equal to or greater than zero if
 * "s1" is lexicographically less than, equal to or greater than "s2".
 */
extern NCBI_XCONNECT_EXPORT
int strncasecmp(const char* s1, const char* s2, size_t n);

#endif/*HAVE_STRCASECMP*/


#ifdef   strupr
#  undef strupr
#  undef strlwr
#endif
#define  strupr        NCBI_strupr
#define  strlwr        NCBI_strlwr

/* Convert a string to uppercase, then return pointer to
 * the altered string. Because the conversion is made in place, the
 * returned pointer is the same as the passed one.
 */
extern NCBI_XCONNECT_EXPORT
char* strupr(char* s);

/* Convert a string to lowercase, then return pointer to
 * the altered string. Because the conversion is made in place, the
 * returned pointer is the same as the passed one.
 */
extern NCBI_XCONNECT_EXPORT
char* strlwr(char* s);


/* Copy not more than "n" characters from string "s2" into "s1"
 * and return the result, which is always null-terminated.
 * NOTE: The difference of this function from standard strncpy() is in
 * that the result is always null-terminated and that the function does not
 * pad "s1" with null bytes should "s2" be shorter than "n" characters.
 */
extern NCBI_XCONNECT_EXPORT
char* strncpy0(char* s1, const char* s2, size_t n);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CONNECT___NCBI_ANSI_EXT__H */
