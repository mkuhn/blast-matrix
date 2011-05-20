#ifndef CONNECT_EXT___NCBI_LOCALNET__H
#define CONNECT_EXT___NCBI_LOCALNET__H

/* $Id: ncbi_localnet.h 186552 2010-03-22 16:11:37Z lavr $
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
 *   Get IP address of a CGI client and determine the IP locality
 *
 *   NOTE:  This is an internal NCBI-only API not for export to third parties.
 *
 */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Return non-zero (true) if called from within a CGI that was invoked
 * by NCBI local client; zero otherwise.
 * 
 * @param
 *     tracking_env - buffer with tracking environment, or if NULL then
 *     the process environment (as provided by the system) will be used.
 */
extern int/*bool*/ NcbiIsLocalCgiClient
(const char* const* tracking_env);


/**
 * Return non-zero (true) if the IP address (in network byte order)
 * provided as an agrument, is a loopback one;  zero otherwise.
 */
extern int/*bool*/ NcbiIsLoopbackIP
(unsigned int ip);


/**
 * Return non-zero (true) if the IP address (in network byte order)
 * provided as an agrument, is a local one;  zero otherwise.
 */
extern int/*bool*/ NcbiIsLocalIP
(unsigned int ip);


/**
 * Return IP address (in network by order) of the CGI client.
 * Result is undetermined if called not from within CGI executable.
 * Return 0 if the IP address cannot be obtained.
 *
 * @param
 *     buf - buffer where client name will be saved, could be NULL.
 * @param
 *     buf_size - buffer size.
 * @param
 *     tracking_env - buffer with tracking environment, or if NULL then
 *     the process environment (as provided by the system) will be used.
 */

typedef enum {
    eCgiClientIP_TryMost  = 0, /* Try most of known environment variables   */
    eCgiClientIP_TryAll   = 1, /* Try all (NI_CLIENT_IPADDR incl.) env. vars*/
    eCgiClientIP_TryLeast = 2  /* Try to detect caller's IP only, not origin*/
} ECgiClientIP;

extern unsigned int NcbiGetCgiClientIPEx
(ECgiClientIP       flag,
 char*              buf,
 size_t             buf_size,
 const char* const* tracking_env
 );

/* NcbiGetCgiClientIPEx(., NULL, 0, .) */
extern unsigned int NcbiGetCgiClientIP
(ECgiClientIP       flag,
 const char* const* tracking_env
 );


#ifdef __cplusplus
}  /* extern "C" */
#endif


#endif /*CONNECT_EXT___NCBI_LOCALNET__H*/
