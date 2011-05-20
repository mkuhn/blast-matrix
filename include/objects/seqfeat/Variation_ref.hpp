/* $Id: Variation_ref.hpp 206953 2010-09-30 16:51:18Z dicuccio $
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

/// @file Variation_ref.hpp
/// User-defined methods of the data storage class.
///
/// This file was originally generated by application DATATOOL
/// using the following specifications:
/// 'seqfeat.asn'.
///
/// New methods or data members can be added to it if needed.
/// See also: Variation_ref_.hpp


#ifndef OBJECTS_SEQFEAT_VARIATION_REF_HPP
#define OBJECTS_SEQFEAT_VARIATION_REF_HPP


// generated includes
#include <objects/seqfeat/Variation_ref_.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

/////////////////////////////////////////////////////////////////////////////
class NCBI_SEQFEAT_EXPORT CVariation_ref : public CVariation_ref_Base
{
    typedef CVariation_ref_Base Tparent;
public:
    // constructor
    CVariation_ref(void);
    // destructor
    ~CVariation_ref(void);

    /// Enum governing sequence types for supplied sequence strings
    enum ESeqType {
        eSeqType_na,
        eSeqType_aa
    };

    /// Set a standard single nucleotide variant.  The replaces set can include
    /// empty strings and/or '-' as a character to indicate a deletion.
    void SetSNV(const vector<string>& replaces,
                ESeqType seq_type);
    bool IsSNV() const;

    /// Set a standard multinucleotide variant.  The replaces set can include
    /// empty strings and/or '-' as a character to indicate a deletion.
    void SetMNP(const vector<string>& replaces,
                ESeqType seq_type);
    bool IsMNP() const;

    /// Make this variant a deletion
    void SetDeletion();
    bool IsDeletion() const;

    /// Make this variant an insertion
    void SetInsertion(const string& sequence, ESeqType seq_type);
    bool IsInsertion() const;

    /// Make this variant an insertion of unknown sequence
    void SetInsertion();

    /// Make this variant an insertion
    void SetDeletionInsertion(const string& sequence, ESeqType seq_type);
    bool IsDeletionInsertion() const;

    /// Set the standard fields for a microsatellite.  This API establishes a
    /// microsatellite with a range of possible observed repeats.
    void SetMicrosatellite(const string& nucleotide_seq,
                           TSeqPos min_repeats, TSeqPos max_repeats);
    bool IsMicrosatellite() const;

    /// Set the standard fields for a microsatellite.  This API establishes a
    /// microsatellite with a fixed set of possible observed repeats
    void SetMicrosatellite(const string& nucleotide_seq,
                           const vector<TSeqPos>& observed_repeats);

    /// Make this variant a copy number variant.  NOTE: This API variant
    /// establishes a CNV of unknown copy number
    void SetCNV();
    bool IsCNV() const;

    /// Special subtype of CNV: 'gain' - an unspecified increase in copy number
    void SetGain();
    bool IsGain() const;

    /// Special subtype of CNV: 'loss' - an unspecified decrease in copy number
    void SetLoss();
    bool IsLoss() const;

    /// Make this variant a copy number variant.  NOTE: This API variant
    /// establishes a CNV with a range of possible copies
    void SetCNV(TSeqPos min_copy, TSeqPos max_copy);

    /// Make this variant a copy number variant.  NOTE: This API variant
    /// establishes a CNV with a fixed set of possible copies
    void SetCNV(const vector<TSeqPos>& observed_copies);

    /// The feature represents an inversion at the specified location
    /// The provided location should be upstream and on the opposite strand
    void SetInversion(const CSeq_loc& other_loc);
    bool IsInversion() const;

    /// The feature represents an eversion at the specified location
    /// The provided location should be downstream and on the opposite strand
    void SetEversion(const CSeq_loc& other_loc);
    bool IsEversion() const;

    /// The feature represents a translocation event
    /// The provided location can be anywhere; a special case exists when the
    /// provided location is on a different chromosome, in which case the
    /// feature is considered a transchromosomal rearrangement
    void SetTranslocation(const CSeq_loc& other_loc);
    bool IsTranslocation() const;

    /// Establish a uniparental disomy mark-up
    void SetUniparentalDisomy();
    bool IsUniparentalDisomy() const;

    /// Create a complex undescribed variant
    void SetComplex();
    bool IsComplex() const;

    /// Create a variant of unknown type
    void SetUnknown();
    bool IsUnknown() const;

    /// Create a variant of type 'other'
    void SetOther();
    bool IsOther() const;

    /// Validate that all semantic fields are correct
    void Validate();

private:
    // Prohibit copy constructor and assignment operator
    CVariation_ref(const CVariation_ref& value);
    CVariation_ref& operator=(const CVariation_ref& value);

};

/////////////////// CVariation_ref inline methods

// constructor
inline
CVariation_ref::CVariation_ref(void)
{
}


/////////////////// end of CVariation_ref inline methods


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE


#endif // OBJECTS_SEQFEAT_VARIATION_REF_HPP
/* Original file checksum: lines: 86, chars: 2492, CRC32: ef8e854f */
