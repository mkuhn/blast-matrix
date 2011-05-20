/* $Id: cleanup.cpp 196740 2010-07-08 13:40:20Z bollin $
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
 * Author:  Robert Smith
 *
 * File Description:
 *   Basic Cleanup of CSeq_entries.
 *
 */
#include <ncbi_pch.hpp>
#include <corelib/ncbistd.hpp>
#include <serial/serialbase.hpp>
#include <objects/seq/Bioseq.hpp>
#include <objects/seq/Seq_annot.hpp>
#include <objects/seqset/Seq_entry.hpp>
#include <objects/seqset/Bioseq_set.hpp>
#include <objects/seqfeat/Seq_feat.hpp>
#include <objects/submit/Seq_submit.hpp>


#include <objmgr/util/sequence.hpp>
#include <objtools/cleanup/cleanup.hpp>

#include "cleanupp.hpp"


BEGIN_NCBI_SCOPE
BEGIN_SCOPE(objects)

enum EChangeType {
    eChange_UNKNOWN
};

// *********************** CCleanup implementation **********************


CCleanup::CCleanup(CScope* scope) : m_Scope( scope )
{
}


CCleanup::~CCleanup(void)
{
}


void CCleanup::SetScope(CScope* scope)
{
    m_Scope = scope;
}


static
CRef<CCleanupChange> makeCleanupChange(Uint4 options)
{
    CRef<CCleanupChange> changes;
    if (! (options  &  CCleanup::eClean_NoReporting)) {
        changes.Reset(new CCleanupChange);        
    }
    return changes;
}

CConstRef<CCleanupChange> CCleanup::BasicCleanup(CSeq_entry& se, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(se);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CSeq_submit& ss, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(ss);
    return changes;
}


/// Cleanup a Bioseq. 
CConstRef<CCleanupChange> CCleanup::BasicCleanup(CBioseq& bs, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(bs);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CBioseq_set& bss, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(bss);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CSeq_annot& sa, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(sa);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CSeq_feat& sf, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(sf);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CSeq_entry_Handle& seh, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(seh);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CBioseq_Handle& bsh,    Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(bsh);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CBioseq_set_Handle& bssh, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(bssh);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CSeq_annot_Handle& sah, Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(sah);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::BasicCleanup(CSeq_feat_Handle& sfh,  Uint4 options)
{
    CRef<CCleanupChange> changes(makeCleanupChange(options));
    CCleanup_imp clean_i(changes, m_Scope, options);
    clean_i.BasicCleanup(sfh);
    return changes;
}





// *********************** Extended Cleanup implementation ********************
CConstRef<CCleanupChange> CCleanup::ExtendedCleanup(const CSeq_entry& se)
{
    CRef<CCleanupChange> changes(makeCleanupChange(0));
    CCleanup_imp clean_i(changes, m_Scope, 0);
    clean_i.ExtendedCleanup(m_Scope->GetSeq_entryHandle(se));
    
    return changes;
}


CConstRef<CCleanupChange> CCleanup::ExtendedCleanup(const CSeq_submit& ss)
{
    CRef<CCleanupChange> changes(makeCleanupChange(0));
    CCleanup_imp clean_i(changes, m_Scope, 0);
    clean_i.ExtendedCleanup(ss);
    return changes;
}


CConstRef<CCleanupChange> CCleanup::ExtendedCleanup(const CSeq_annot& sa)
{
    CRef<CCleanupChange> changes(makeCleanupChange(0));
    CCleanup_imp clean_i(changes, m_Scope, 0);
    clean_i.ExtendedCleanup(m_Scope->GetSeq_annotHandle(sa));
    return changes;
}


// *********************** CCleanupChange implementation **********************


CCleanupChange::CCleanupChange()
{
}


size_t CCleanupChange::ChangeCount() const
{
    return m_Changes.count();
}


bool CCleanupChange::IsChanged(CCleanupChange::EChanges e) const
{
    return m_Changes.test(e);
}


void CCleanupChange::SetChanged(CCleanupChange::EChanges e)
{
    m_Changes.set(e);
}


vector<CCleanupChange::EChanges> CCleanupChange::GetAllChanges() const
{
    vector<EChanges>  result;
    for (size_t i = eNoChange + 1; i < m_Changes.size(); ++i) {
        if (m_Changes.test(i)) {
            result.push_back( (EChanges) i);
        }
    }
    return result;
}


vector<string> CCleanupChange::GetAllDescriptions() const
{
    vector<string>  result;
    for (size_t i = eNoChange + 1; i < m_Changes.size(); ++i) {
        if (m_Changes.test(i)) {
            result.push_back( GetDescription((EChanges) i) );
        }
    }
    return result;
}


string CCleanupChange::GetDescription(EChanges e)
{
    if (e <= eNoChange  ||  e >= eNumberofChangeTypes) {
        return sm_ChangeDesc[eNoChange];
    }
    return sm_ChangeDesc[e];
}


const char* CCleanupChange::sm_ChangeDesc[] = {
    "Invalid Change Code",
    // set when strings are changed.
    "Trim Spaces",
    "Clean Double Quotes",
    // set when lists are sorted or uniqued.
    "Clean Qualifiers List",
    "Clean Dbxrefs List",
    "Clean CitonFeat List",
    "Clean Keywords List",
    "Clean Subsource List",
    "Clean Orgmod List",
    // Set when fields are moved or have content changes
    "Repair BioseqMol",
    "Change Feature Key", //10
    "Normalize Authors",
    "Change Publication",
    "Change Qualifiers",
    "Change Dbxrefs",
    "Change Keywords",
    "Change Subsource",
    "Change Orgmod",
    "Change Exception",
    "Change Comment",
    // Set when fields are rescued
    "Change tRna", //20
    "Change rRna",
    "Change ITS",
    "Change Anticodon",
    "Change Code Break",
    "Change Genetic Code",
    "Copy GeneXref",
    "Copy ProtXref",
    // set when locations are repaired
    "Change Seqloc",
    "Change Strand",
    "Change WholeLocation", //30
    // set when MolInfo descriptors are affected
    "Change MolInfo Descriptor",
    // set when prot-xref is removed
    "Remove ProtXref",
    // set when gene-xref is removed
    "Remove GeneXref",
    // set when protein feature is added
    "Add Protein Feature",
    // set when feature is removed
    "Remove Feature",
    // set when feature is moved
    "Move Feature",
    // set when qualifier is removed
    "Remove Qualifier",
    // set when Gene Xref is created
    "Add GeneXref",
    // set when descriptor is removed
    "Remove Descriptor",
    "Remove Keyword", //40
    "Add Descriptor",
    "Move Descriptor",
    "Convert Feature to Descriptor", 
    "Collapse Set",
    "Change Feature Location",
    "Remove Annotation",
    "Convert Feature",
    "Remove Comment",
    "Add BioSource OrgMod",
    "Add BioSource SubSource", //50
    "Change BioSource Genome", 
    "Change BioSource Origin", 
    "Change BioSource Other",
    "Change SeqId", 
    "Remove Empty Publication",
    "Add Qualifier",
    "Cleanup Date",
    "Change BioseqInst",
    "Remove SeqID",
    "Add ProtXref",
    "Change Partial",
    "Change Prot Names",
    // set when any other change is made.
    "Change Other", 
    "Invalid Change Code"
};

END_SCOPE(objects)
END_NCBI_SCOPE
