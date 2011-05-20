/* $Id: Dense_diag.hpp 207227 2010-10-04 15:26:24Z grichenk $
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
 */

/// @file Dense_diag.hpp
/// User-defined methods of the data storage class.
///
/// This file was originally generated by application DATATOOL
/// using the following specifications:
/// 'seqalign.asn'.
///
/// New methods or data members can be added to it if needed.
/// See also: Dense_diag_.hpp


#ifndef OBJECTS_SEQALIGN_DENSE_DIAG_HPP
#define OBJECTS_SEQALIGN_DENSE_DIAG_HPP


// generated includes
#include <objects/seqalign/Dense_diag_.hpp>

#include <objects/seqalign/seqalign_exception.hpp>
#include <util/range.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

class CSeq_interval;

/////////////////////////////////////////////////////////////////////////////
class NCBI_SEQALIGN_EXPORT CDense_diag : public CDense_diag_Base
{
    typedef CDense_diag_Base Tparent;
public:
    // constructor
    CDense_diag(void);
    // destructor
    ~CDense_diag(void);

    /// Validators
    TDim    CheckNumRows(void) const;
    void    Validate    ()     const;

    /// GetSeqRange
    CRange<TSeqPos> GetSeqRange (TDim row) const;
    TSeqPos         GetSeqStart (TDim row) const;
    TSeqPos         GetSeqStop  (TDim row) const;
    ENa_strand      GetSeqStrand(TDim row) const;

    /// Offset row's coords
    void OffsetRow(TDim row, TSignedSeqPos offset);

    CRef<CSeq_interval> CreateRowSeq_interval(TDim row) const;

private:
    // Prohibit copy constructor and assignment operator
    CDense_diag(const CDense_diag& value);
    CDense_diag& operator=(const CDense_diag& value);

};

/////////////////// CDense_diag inline methods

// constructor
inline
CDense_diag::CDense_diag(void)
{
}


inline
TSeqPos CDense_diag::GetSeqStart(TDim row) const
{
    if (row < 0  ||  row >= GetDim()) {
        NCBI_THROW(CSeqalignException, eInvalidRowNumber,
                   "CDense_diag::GetSeqStart():"
                   " Invalid row number");
    }
    return GetStarts()[row];
}


inline
TSeqPos CDense_diag::GetSeqStop(TDim row) const
{
    return GetSeqStart(row) + GetLen();
}


inline
CRange<TSeqPos> CDense_diag::GetSeqRange(TDim row) const
{
    return CRange<TSeqPos>(GetSeqStart(row), GetSeqStop(row));
}


inline
CDense_diag::TDim CDense_diag::CheckNumRows() const
{
    const size_t& dim = GetDim();
    if (dim != GetIds().size()  ||  dim != GetStarts().size()) {
        NCBI_THROW(CSeqalignException, eInvalidAlignment,
                   "CDense_diag::CheckNumRows()"
                   " dim is not consistent with ids.size & starts.size");
    }
    return dim;
}


/////////////////// end of CDense_diag inline methods


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE

#endif // OBJECTS_SEQALIGN_DENSE_DIAG_HPP
/* Original file checksum: lines: 94, chars: 2608, CRC32: ba2ec5bc */
