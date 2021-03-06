/* $Id: Seq_descr.cpp 193581 2010-06-04 17:24:24Z gouriano $
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
 * Author:  .......
 *
 * File Description:
 *   .......
 *
 * Remark:
 *   This code was originally generated by application DATATOOL
 *   using the following specifications:
 *   'seq.asn'.
 */

// standard includes
#include <ncbi_pch.hpp>
#include <objects/seq/Seqdesc.hpp>

// generated includes
#include <objects/seq/Seq_descr.hpp>

#include <corelib/ncbi_param.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

// destructor
CSeq_descr::~CSeq_descr(void)
{
}

NCBI_PARAM_DECL(bool, OBJECTS, SEQ_DESCR_ALLOW_EMPTY);
NCBI_PARAM_DEF_EX(bool, OBJECTS, SEQ_DESCR_ALLOW_EMPTY, false,
                  eParam_NoThread, OBJECTS_SEQ_DESCR_ALLOW_EMPTY);

void CSeq_descr::PostRead(void) const
{
    static NCBI_PARAM_TYPE(OBJECTS, SEQ_DESCR_ALLOW_EMPTY) sx_Value;
    if ( !sx_Value.Get() && Get().empty() ) {
        NCBI_THROW(CSerialException, eInvalidData,
                   "empty Seq-descr is not allowed");
    }
}

void CSeq_descr::PreWrite(void) const
{
    static NCBI_PARAM_TYPE(OBJECTS, SEQ_DESCR_ALLOW_EMPTY) sx_Value;
    if ( !sx_Value.Get() && Get().empty() ) {
        NCBI_THROW(CSerialException, eInvalidData,
                   "empty Seq-descr is not allowed");
    }
}

END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE

/* Original file checksum: lines: 65, chars: 1883, CRC32: ed475d39 */
