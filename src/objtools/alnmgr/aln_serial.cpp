/*  $Id: aln_serial.cpp 162058 2009-06-02 20:06:14Z todorov $
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
* Author:  Kamen Todorov, NCBI
*
* File Description:
*   Alignments Serialization
*
* ===========================================================================
*/


#include <ncbi_pch.hpp>
#include <objtools/alnmgr/aln_serial.hpp>


BEGIN_NCBI_SCOPE
USING_SCOPE(objects);


ostream& operator<<(ostream& out, const CPairwiseAln::TRng& rng)
{
    if (rng.GetFrom() < rng.GetToOpen()) {
        out << "["
            << rng.GetFrom() << ", " 
            << rng.GetToOpen()
            << ")"
            << " len: "
            << rng.GetLength();
    } else {
        out << "<"
            << rng.GetFrom() << ", " 
            << rng.GetTo()
            << ">"
            << " len: "
            << rng.GetLength();
    }
    return out;
}


ostream& operator<<(ostream& out, const IAlnSegment::ESegTypeFlags& flags)
{
    return out
        << (flags & IAlnSegment::fAligned ? "fAligned " : "")
        << (flags & IAlnSegment::fGap ? "fGap " : "")
        << (flags & IAlnSegment::fReversed ? "fReversed " : "")
        << (flags & IAlnSegment::fInvalid ? "fInvalid " : "");
}


ostream& operator<<(ostream& out, const IAlnSegment& aln_seg)
{
    return out
        << " Anchor Rng: " << aln_seg.GetAlnRange()
        << " Rng: " << aln_seg.GetRange()
        << " type: " << (IAlnSegment::ESegTypeFlags) aln_seg.GetType();
}


ostream& operator<<(ostream& out, const CPairwiseAln::TAlnRng& aln_rng)
{
    return out 
        << "["
        << aln_rng.GetFirstFrom() << ", " 
        << aln_rng.GetSecondFrom() << ", "
        << aln_rng.GetLength() << ", " 
        << (aln_rng.IsDirect() ? "direct" : "reverse") 
        << "]";
}


ostream& operator<<(ostream& out, const CPairwiseAln::EFlags& flags)
{
    out << " Flags = " << NStr::UIntToString(flags, 0, 2)
        << ":" << endl;
    
    if (flags & CPairwiseAln::fKeepNormalized) out << "fKeepNormalized" << endl;
    if (flags & CPairwiseAln::fAllowMixedDir) out << "fAllowMixedDir" << endl;
    if (flags & CPairwiseAln::fAllowOverlap) out << "fAllowOverlap" << endl;
    if (flags & CPairwiseAln::fAllowAbutting) out << "fAllowAbutting" << endl;
    if (flags & CPairwiseAln::fNotValidated) out << "fNotValidated" << endl;
    if (flags & CPairwiseAln::fInvalid) out << "fInvalid" << endl;
    if (flags & CPairwiseAln::fUnsorted) out << "fUnsorted" << endl;
    if (flags & CPairwiseAln::fDirect) out << "fDirect" << endl;
    if (flags & CPairwiseAln::fReversed) out << "fReversed" << endl;
    if ((flags & CPairwiseAln::fMixedDir) == CPairwiseAln::fMixedDir) out << "fMixedDir" << endl;
    if (flags & CPairwiseAln::fOverlap) out << "fOverlap" << endl;
    if (flags & CPairwiseAln::fAbutting) out << "fAbutting" << endl;

    return out;
}


ostream& operator<<(ostream& out, const TAlnSeqIdIRef& aln_seq_id_iref)
{
    out << aln_seq_id_iref->AsString()
        << " (base_width=" << aln_seq_id_iref->GetBaseWidth() 
        << ")";
    return out;
}        


ostream& operator<<(ostream& out, const CPairwiseAln& pairwise_aln)
{
    out << "CPairwiseAln between ";

    out << pairwise_aln.GetFirstId() << " and "
        << pairwise_aln.GetSecondId();
    
    cout << " with flags=" << pairwise_aln.GetFlags() << " and segments:" << endl;

    ITERATE (CPairwiseAln, aln_rng_it, pairwise_aln) {
        out << *aln_rng_it;
    }
    return out << endl;
}


ostream& operator<<(ostream& out, const CAnchoredAln& anchored_aln)
{
    out << "CAnchorAln has score of " << anchored_aln.GetScore() << " and contains " 
        << anchored_aln.GetDim() << " pair(s) of rows:" << endl;
    ITERATE(CAnchoredAln::TPairwiseAlnVector, pairwise_aln_i, anchored_aln.GetPairwiseAlns()) {
        out << **pairwise_aln_i;
    }
    return out << endl;
}


ostream& operator<<(ostream& out, const TAlnStats& aln_stats)
{
    out << "Number of alignments: " << aln_stats.GetAlnCount() << endl;
    out << "IsCanonicalQueryAnchored: " << aln_stats.IsCanonicalQueryAnchored() << endl;
    out << "IsCanonicalMultiple: " << aln_stats.IsCanonicalMultiple() << endl;
    out << "CanBeAnchored: " << aln_stats.CanBeAnchored() << endl;
    out << endl;
    out << "IdVec (" << aln_stats.GetIdVec().size() << "):" << endl;
    ITERATE(TAlnStats::TIdVec, it, aln_stats.GetIdVec()) {
        out << (*it)->AsString() << " (base_width=" << (*it)->GetBaseWidth() << ")" << endl;
    }
    out << endl;
    out << "IdMap (" << aln_stats.GetIdMap().size() << "):" << endl;
    ITERATE(TAlnStats::TIdMap, it, aln_stats.GetIdMap()) {
        out << it->first->AsString() << " (base_width=" << it->first->GetBaseWidth() << ")" << endl;
    }
    out << endl;
    out << "AnchorIdVec (" << aln_stats.GetAnchorIdVec().size() << "):" << endl;
    ITERATE(TAlnStats::TIdVec, it, aln_stats.GetAnchorIdVec()) {
        out << (*it)->AsString() << " (base_width=" << (*it)->GetBaseWidth() << ")" << endl;
    }
    out << endl;
    out << "AnchorIdMap (" << aln_stats.GetAnchorIdMap().size() << "):" << endl;
    ITERATE(TAlnStats::TIdMap, it, aln_stats.GetAnchorIdMap()) {
        out << it->first->AsString() << " (base_width=" << it->first->GetBaseWidth() << ")" << endl;
    }
    out << endl;
    out << "AnchorIdxVec (" << aln_stats.GetAnchorIdxVec().size() << "):" << endl;
    ITERATE(TAlnStats::TIdxVec, it, aln_stats.GetAnchorIdxVec()) {
        out << *it << endl;
    }
    out << endl;
    for (size_t aln_idx = 0;  aln_idx < aln_stats.GetAlnCount();  ++aln_idx) {
        TAlnStats::TDim dim = aln_stats.GetDimForAln(aln_idx);
        out << "Alignment " << aln_idx << " has " 
            << dim << " rows:" << endl;
        for (TAlnStats::TDim row = 0;  row < dim;  ++row) {
            out << aln_stats.GetSeqIdsForAln(aln_idx)[row]->AsString();
            out << endl;
        }
        out << endl;
    }
    return out;
}


END_NCBI_SCOPE


