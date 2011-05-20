/* $Id: Variation_ref.cpp 206953 2010-09-30 16:51:18Z dicuccio $
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
 *   'seqfeat.asn'.
 */

// standard includes
#include <ncbi_pch.hpp>

// generated includes
#include <objects/seqfeat/Variation_ref.hpp>
#include <objects/seq/Seq_literal.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

// destructor
CVariation_ref::~CVariation_ref(void)
{
}

static void s_SetReplaces(CVariation_ref& ref,
                          const vector<string>& replaces,
                          CVariation_ref::ESeqType seq_type,
                          CVariation_inst::EType var_type)
{
    list< CRef<CDelta_item> > items;
    bool has_del = false;

    ITERATE (vector<string>, it, replaces) {
        string rep(*it);
        NStr::ToUpper(rep);
        NStr::TruncateSpacesInPlace(rep);

        if (rep.empty()  ||  rep == "-") {
            has_del = true;
        } else {
            CRef<CDelta_item> item(new CDelta_item);
            item->SetSeq().SetLiteral().SetLength(rep.size());
            if (seq_type == CVariation_ref::eSeqType_na) {
                item->SetSeq().SetLiteral().SetSeq_data().SetIupacna().Set(rep);
            } else {
                item->SetSeq().SetLiteral().SetSeq_data().SetIupacaa().Set(rep);
            }
            items.push_back(item);
        }
    }

    if (has_del  &&  items.size()) {
        ///
        /// both deletion and replaces
        /// therefore, we are a complex set
        ///
        CRef<CVariation_ref> sub;

        ref.SetData().SetSet().SetType
            (CVariation_ref::TData::TSet::eData_set_type_compound);
        
        /// deletion first
        sub.Reset(new CVariation_ref);
        sub->SetData().SetInstance().SetType(CVariation_inst::eType_del);
        sub->SetData().SetInstance().SetDelta().clear();
        ref.SetData().SetSet().SetVariations().push_back(sub);

        /// then the replaces
        sub.Reset(new CVariation_ref);
        sub->SetData().SetInstance().SetType(var_type);
        sub->SetData().SetInstance().SetDelta()
            .insert(sub->SetData().SetInstance().SetDelta().end(),
                    items.begin(), items.end());
        ref.SetData().SetSet().SetVariations().push_back(sub);
    }
    else if (has_del) {
        ref.SetData().SetInstance().SetDelta().clear();
    }
    else if (items.size()) {
        ref.SetData().SetInstance().SetType(var_type);
        ref.SetData().SetInstance().SetDelta()
            .insert(ref.SetData().SetInstance().SetDelta().end(),
                    items.begin(), items.end());
    }

    /**
    ITERATE (vector<string>, it, replaces) {
        string rep(*it);
        NStr::ToUpper(rep);

        CRef<CVariation_ref> ref(new CVariation_ref);
        CVariation_inst& inst = ref->SetData().SetInstance();
        inst.SetType(CVariation_inst::eType_snp);
        inst.SetDelta().clear();

        CRef<CDelta_item> item(new CDelta_item);
        item->SetSeq().SetLiteral().SetLength(rep.size());
        if (seq_type == eSeqType_na) {
            item->SetSeq().SetLiteral().SetSeq_data().SetIupacna().Set(rep);
        } else {
            item->SetSeq().SetLiteral().SetSeq_data().SetIupacaa().Set(rep);
        }
        inst.SetDelta().push_back(item);

        SetData().SetSet().SetVariations().push_back(ref);
    }
    SetData().SetSet().SetType(CVariation_ref::TData::TSet::eData_set_type_population);
    **/
}


void CVariation_ref::SetSNV(const vector<string>& replaces,
                            ESeqType seq_type)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetDelta().clear();

    s_SetReplaces(*this, replaces, seq_type,
                  CVariation_inst::eType_snv);
}


bool CVariation_ref::IsSNV() const
{
    if (GetData().IsInstance()  &&
        GetData().GetInstance().IsSetType()  &&
        GetData().GetInstance().GetType() == CVariation_inst::eType_snv) {
        return true;
    }
    if (GetData().IsSet()) {
        ITERATE (TData::TSet::TVariations, it, GetData().GetSet().GetVariations()) {
            const CVariation_ref& ref = **it;
            if (ref.GetData().IsInstance()  &&
                ref.GetData().GetInstance().IsSetType()  &&
                ref.GetData().GetInstance().GetType() == CVariation_inst::eType_snv) {
                return true;
            }
        }
    }

    return false;
}


void CVariation_ref::SetMNP(const vector<string>& replaces,
                            ESeqType seq_type)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetDelta().clear();

    s_SetReplaces(*this, replaces, seq_type,
                  CVariation_inst::eType_mnp);
}

bool CVariation_ref::IsMNP() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_mnp;
}



void CVariation_ref::SetDeletion()
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetDelta().clear();

    inst.SetType(CVariation_inst::eType_del);
}

bool CVariation_ref::IsDeletion() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_del;
}



void CVariation_ref::SetInsertion(const string& sequence, ESeqType seq_type)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetDelta().clear();

    CRef<CDelta_item> item(new CDelta_item);
    item->SetSeq().SetThis();
    inst.SetDelta().push_back(item);

    vector<string> replaces;
    replaces.push_back(sequence);
    s_SetReplaces(*this, replaces, seq_type,
                  CVariation_inst::eType_ins);
}

bool CVariation_ref::IsInsertion() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_ins;
}



void CVariation_ref::SetInsertion()
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_ins);

    CRef<CDelta_item> item(new CDelta_item);
    item->SetAction(CDelta_item::eAction_ins_before);
    inst.SetDelta().clear();
    inst.SetDelta().push_back(item);
}


void CVariation_ref::SetDeletionInsertion(const string& sequence,
                                          ESeqType seq_type)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetDelta().clear();

    CRef<CDelta_item> item;

    item.Reset(new CDelta_item);
    item->SetAction(CDelta_item::eAction_del_at);
    inst.SetDelta().push_back(item);

    vector<string> replaces;
    replaces.push_back(sequence);
    s_SetReplaces(*this, replaces, seq_type,
                  CVariation_inst::eType_delins);
}

bool CVariation_ref::IsDeletionInsertion() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_delins;
}



void CVariation_ref::SetMicrosatellite(const string& nucleotide_seq,
                                       TSeqPos min_repeats,
                                       TSeqPos max_repeats)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetDelta().clear();

    vector<string> replaces;
    replaces.push_back(nucleotide_seq);
    s_SetReplaces(*this, replaces, eSeqType_na,
                  CVariation_inst::eType_microsatellite);

    inst.SetDelta().front()->SetMultiplier(min_repeats);
    inst.SetDelta().front()
        ->SetMultiplier_fuzz().SetRange().SetMin(min_repeats);
    inst.SetDelta().front()
        ->SetMultiplier_fuzz().SetRange().SetMax(max_repeats);
}

bool CVariation_ref::IsMicrosatellite() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_microsatellite;
}



void CVariation_ref::SetMicrosatellite(const string& nucleotide_seq,
                                       const vector<TSeqPos>& observed_repeats)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetDelta().clear();

    vector<string> replaces;
    replaces.push_back(nucleotide_seq);
    s_SetReplaces(*this, replaces, eSeqType_na,
                  CVariation_inst::eType_microsatellite);

    inst.SetDelta().front()->SetMultiplier(observed_repeats.front());
    if (observed_repeats.size() > 1) {
        std::copy(observed_repeats.begin(),
                  observed_repeats.end(),
                  back_inserter(inst.SetDelta().front()
                                ->SetMultiplier_fuzz().SetAlt()));
    }
}


void CVariation_ref::SetCNV()
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_cnv);
    inst.SetDelta().clear();

    CRef<CDelta_item> item(new CDelta_item);
    item->SetSeq().SetThis();
    item->SetMultiplier_fuzz().SetLim(CInt_fuzz::eLim_unk);

    inst.SetDelta().push_back(item);
}

bool CVariation_ref::IsCNV() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_cnv;
}



void CVariation_ref::SetGain()
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_cnv);
    inst.SetDelta().clear();

    CRef<CDelta_item> item(new CDelta_item);
    item->SetSeq().SetThis();
    item->SetMultiplier_fuzz().SetLim(CInt_fuzz::eLim_gt);

    inst.SetDelta().push_back(item);
}

bool CVariation_ref::IsGain() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_cnv  &&
           GetData().GetInstance().IsSetDelta()  &&
           GetData().GetInstance().GetDelta().size()  &&
           GetData().GetInstance().GetDelta().front()->IsSetMultiplier_fuzz()  &&
           GetData().GetInstance().GetDelta().front()->GetMultiplier_fuzz().IsLim()  &&
           GetData().GetInstance().GetDelta().front()->GetMultiplier_fuzz().GetLim() == CInt_fuzz::eLim_gt;

}



void CVariation_ref::SetLoss()
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_cnv);
    inst.SetDelta().clear();

    CRef<CDelta_item> item(new CDelta_item);
    item->SetSeq().SetThis();
    item->SetMultiplier_fuzz().SetLim(CInt_fuzz::eLim_lt);

    inst.SetDelta().push_back(item);
}

bool CVariation_ref::IsLoss() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_cnv  &&
           GetData().GetInstance().IsSetDelta()  &&
           GetData().GetInstance().GetDelta().size()  &&
           GetData().GetInstance().GetDelta().front()->IsSetMultiplier_fuzz()  &&
           GetData().GetInstance().GetDelta().front()->GetMultiplier_fuzz().IsLim()  &&
           GetData().GetInstance().GetDelta().front()->GetMultiplier_fuzz().GetLim() == CInt_fuzz::eLim_lt;

}


void CVariation_ref::SetCNV(TSeqPos min_copy, TSeqPos max_copy)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_cnv);
    inst.SetDelta().clear();

    CRef<CDelta_item> item(new CDelta_item);
    item->SetSeq().SetThis();
    item->SetMultiplier_fuzz().SetRange().SetMin(min_copy);
    item->SetMultiplier_fuzz().SetRange().SetMax(max_copy);

    inst.SetDelta().push_back(item);
}


void CVariation_ref::SetCNV(const vector<TSeqPos>& observed_copies)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_cnv);
    inst.SetDelta().clear();

    CRef<CDelta_item> item(new CDelta_item);
    item->SetSeq().SetThis();
    std::copy(observed_copies.begin(), observed_copies.end(),
              back_inserter(item->SetMultiplier_fuzz().SetAlt()));

    inst.SetDelta().push_back(item);
}

void CVariation_ref::SetInversion(const CSeq_loc& other_loc)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_inverted_copy);
    inst.SetDelta().clear();

    CRef<CDelta_item> item(new CDelta_item);
    item->SetSeq().SetLoc().Assign(other_loc);
    inst.SetDelta().push_back(item);
}

bool CVariation_ref::IsInversion() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_inverted_copy;
}



void CVariation_ref::SetEversion(const CSeq_loc& other_loc)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_everted_copy);
    inst.SetDelta().clear();

    CRef<CDelta_item> item(new CDelta_item);
    item->SetSeq().SetLoc().Assign(other_loc);
    inst.SetDelta().push_back(item);
}

bool CVariation_ref::IsEversion() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_everted_copy;
}



/// The feature represents an eversion at the specified location
/// The provided location can be anywhere; a special case exists when the
/// provided location is on a different chromosome, in which case the
/// feature is considered a transchromosomal rearrangement
void CVariation_ref::SetTranslocation(const CSeq_loc& other_loc)
{
    CVariation_inst& inst = SetData().SetInstance();
    inst.SetType(CVariation_inst::eType_translocation);
    inst.SetDelta().clear();

    CRef<CDelta_item> item;
    item.Reset(new CDelta_item);
    item->SetAction(CDelta_item::eAction_del_at);
    inst.SetDelta().push_back(item);

    item.Reset(new CDelta_item);
    item->SetSeq().SetLoc().Assign(other_loc);
    inst.SetDelta().push_back(item);

}

bool CVariation_ref::IsTranslocation() const
{
    return GetData().IsInstance()  &&
           GetData().GetInstance().IsSetType()  &&
           GetData().GetInstance().GetType() == CVariation_inst::eType_translocation;
}



/// Establish a uniparental disomy mark-up
void CVariation_ref::SetUniparentalDisomy()
{
    SetData().SetUniparental_disomy();
}

bool CVariation_ref::IsUniparentalDisomy() const
{
    return GetData().IsUniparental_disomy();
}



/// Establish a complex undescribed variant
void CVariation_ref::SetComplex()
{
    SetData().SetComplex();
}

bool CVariation_ref::IsComplex() const
{
    return GetData().IsComplex();
}


void CVariation_ref::SetUnknown()
{
    SetData().SetUnknown();
}

bool CVariation_ref::IsUnknown() const
{
    return GetData().IsUnknown();
}


void CVariation_ref::SetOther()
{
    SetData().SetSet().SetType
        (CVariation_ref::TData::TSet::eData_set_type_other);
    SetData().SetSet().SetVariations();

}

bool CVariation_ref::IsOther() const
{
    return GetData().IsSet()  &&
        GetData().GetSet().GetType() ==
        CVariation_ref::TData::TSet::eData_set_type_other;
}



///
/// perform necessary validation
///
void CVariation_ref::Validate()
{
}


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE

/* Original file checksum: lines: 57, chars: 1736, CRC32: 48a35aa2 */
