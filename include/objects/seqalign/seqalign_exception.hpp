/* $Id: seqalign_exception.hpp 150644 2009-01-28 02:53:46Z todorov $
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
 * Author: Philip Johnson
 *
 * File Description: CSeqalignException class
 *
 */

#ifndef OBJECTS_SEQALIGN_SEQALIGN_EXCEPTION_HPP
#define OBJECTS_SEQALIGN_SEQALIGN_EXCEPTION_HPP

#include <corelib/ncbiexpt.hpp>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(objects) // namespace ncbi::objects::

class CSeqalignException : EXCEPTION_VIRTUAL_BASE public CException
{
public:
    enum EErrCode {
        eUnsupported,
        eInvalidAlignment,
        eInvalidInputAlignment,
        eInvalidRowNumber,
        eOutOfRange,
        eInvalidInputData,
        eInvalidSeqId
    };

    virtual const char* GetErrCodeString(void) const
    {
        switch (GetErrCode()) {
        case eUnsupported:           return "eUnsupported";
        case eInvalidAlignment:      return "eInvalidAlignment";
        case eInvalidInputAlignment: return "eInvalidInputAlignment";
        case eInvalidRowNumber:      return "eInvalidRowNumber";
        case eOutOfRange:            return "eOutOfRange";
        case eInvalidInputData:      return "eInvalidInputData";
        case eInvalidSeqId:          return "eInvalidSeqId";
        default:                     return CException::GetErrCodeString();
        }
    }

    NCBI_EXCEPTION_DEFAULT(CSeqalignException, CException);
};

#define _SEQALIGN_ASSERT(expr) \
    do {                                                               \
        if ( !(expr) ) {                                               \
            _ASSERT(expr);                                             \
            NCBI_THROW(CSeqalignException, eInvalidAlignment,          \
                       string("Assertion failed: ") + #expr);          \
        }                                                              \
    } while ( 0 )

END_SCOPE(objects)
END_NCBI_SCOPE

#endif // OBJECTS_SEQALIGN_SEQALIGN_EXCEPTION_HPP
