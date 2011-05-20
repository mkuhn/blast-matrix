/*  $Id: blastdb_aux.hpp 168572 2009-08-18 14:26:04Z camacho $
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
 * Author: Christiam Camacho
 *
 */

/** @file blastdb_aux.hpp
 * Auxiliary functions for BLAST database application
 */

#ifndef __SRC_APP_BLASTDB__BLASTDB_AUX__HPP__
#define __SRC_APP_BLASTDB__BLASTDB_AUX__HPP__

#include <objtools/blast/seqdb_reader/seqdb.hpp>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

/// Convert a string to a CSeqDB ESeqType object
CSeqDB::ESeqType ParseTypeString(const string& str);

END_SCOPE(blast)
END_NCBI_SCOPE

#endif /* __SRC_APP_BLASTDB__BLASTDB_AUX__HPP__ */

