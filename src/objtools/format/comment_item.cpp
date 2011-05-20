/*  $Id: comment_item.cpp 240010 2011-02-02 16:25:17Z rafanovi $
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
* Author:  Mati Shomrat, NCBI
*
* File Description:
*   flat-file generator -- comment item implementation
*
*/
#include <ncbi_pch.hpp>
#include <corelib/ncbistd.hpp>

#include <objects/seqfeat/Seq_feat.hpp>
#include <objects/seq/Seq_hist.hpp>
#include <objects/seq/Seq_hist_rec.hpp>
#include <objects/seq/Seqdesc.hpp>
#include <objects/seq/MolInfo.hpp>
#include <objects/seqfeat/BioSource.hpp>
#include <objects/seqfeat/Org_ref.hpp>
#include <objects/seqfeat/SubSource.hpp>
#include <objects/seqfeat/OrgName.hpp>
#include <objects/seqfeat/OrgMod.hpp>
#include <objects/general/User_object.hpp>
#include <objects/general/User_field.hpp>
#include <objects/general/Object_id.hpp>
#include <objects/general/Date.hpp>
#include <objects/general/Dbtag.hpp>
#include <objmgr/seqdesc_ci.hpp>

#include <objtools/format/formatter.hpp>
#include <objtools/format/text_ostream.hpp>
#include <objtools/format/items/comment_item.hpp>
#include <objtools/format/context.hpp>
#include "utils.hpp"


BEGIN_NCBI_SCOPE
BEGIN_SCOPE(objects)


// static variables initialization
bool CCommentItem::sm_FirstComment = true;

static const string kRefSeq = "REFSEQ";
static const string kRefSeqLink = "<a href=http://www.ncbi.nlm.nih.gov/RefSeq/>REFSEQ</a>";

/////////////////////////////////////////////////////////////////////////////
//
//  CCommentItem

CCommentItem::CCommentItem(CBioseqContext& ctx, bool need_period) :
    CFlatItem(&ctx), m_First(false), m_NeedPeriod(need_period),
    m_CommentInternalIndent(0)
{
    swap(m_First, sm_FirstComment);
}


CCommentItem::CCommentItem
(const string& comment,
 CBioseqContext& ctx,
 const CSerialObject* obj) :
    CFlatItem(&ctx),
    m_First(false), m_NeedPeriod(true),
    m_CommentInternalIndent(0)
{
    m_Comment.push_back( comment );
    ExpandTildes(m_Comment.back(), eTilde_comment);
    swap(m_First, sm_FirstComment);
    if ( obj != 0 ) {
        x_SetObject(*obj);
    }
}

    
CCommentItem::CCommentItem(const CSeqdesc&  desc, CBioseqContext& ctx) :
    CFlatItem(&ctx), m_First(false), m_NeedPeriod(true),
    m_CommentInternalIndent(0)
{
    swap(m_First, sm_FirstComment);
    x_SetObject(desc);
    x_GatherInfo(ctx);
    if ( x_IsCommentEmpty() ) {
        x_SetSkip();
    }
}


CCommentItem::CCommentItem(const CSeq_feat& feat, CBioseqContext& ctx) :
    CFlatItem(&ctx), m_First(false), m_NeedPeriod(true),
    m_CommentInternalIndent(0)
{
    swap(m_First, sm_FirstComment);
    x_SetObject(feat);
    x_GatherInfo(ctx);
    NON_CONST_ITERATE( list<string>, it, m_Comment ) {
        TrimSpacesAndJunkFromEnds( *it );
    }
    if ( x_IsCommentEmpty() ) {
        x_SetSkip();
    }       
}


void CCommentItem::Format
(IFormatter& formatter,
 IFlatTextOStream& text_os) const
{
    formatter.FormatComment(*this, text_os);
}


void CCommentItem::AddPeriod(void)
{
    if( ! m_Comment.empty() ) {
        ncbi::objects::AddPeriod(m_Comment.back());
    }
}


const string& CCommentItem::GetNsAreGapsStr(void)
{
    static const string kNsAreGaps = "The strings of n's in this record represent " \
        "gaps between contigs, and the length of each string corresponds " \
        "to the length of the gap.";
    return kNsAreGaps;
}


string CCommentItem::GetStringForTPA
(const CUser_object& uo,
 CBioseqContext& ctx)
{
    static const string tpa_string = 
        "THIRD PARTY ANNOTATION DATABASE: This TPA record uses data from DDBJ/EMBL/GenBank ";

    if ( !ctx.IsTPA()  ||  ctx.IsRefSeq() ) {
        return kEmptyStr;
    }
    if ( !uo.CanGetType()  ||  !uo.GetType().IsStr()  ||  
         uo.GetType().GetStr() != "TpaAssembly" ) {
        return kEmptyStr;
    }
    
    CBioseq_Handle& seq = ctx.GetHandle();
    if (seq.IsSetInst_Hist()  &&  seq.GetInst_Hist().IsSetAssembly()) {
        return kEmptyStr;
    }

    string id;
    vector<string> accessions;
    ITERATE (CUser_object::TData, curr, uo.GetData()) {
        const CUser_field& uf = **curr;
        if ( !uf.CanGetData()  ||  !uf.GetData().IsFields() ) {
            continue;
        }

        ITERATE (CUser_field::C_Data::TFields, ufi, uf.GetData().GetFields()) {
            if( !(*ufi)->CanGetData()  ||  !(*ufi)->GetData().IsStr()  ||
                !(*ufi)->CanGetLabel() ) {
                continue;
            }
            const CObject_id& oid = (*ufi)->GetLabel();
            if ( oid.IsStr()  &&  
                 (NStr::CompareNocase(oid.GetStr(), "accession") == 0) ) {
                string acc = (*ufi)->GetData().GetStr();
                if ( !acc.empty() ) {
                    accessions.push_back(NStr::ToUpper(acc));
                }
            }
        }
    }
    if ( accessions.empty() ) {
        return kEmptyStr;
    }

    CNcbiOstrstream text;
    text << tpa_string << ((accessions.size() > 1) ? "entries " : "entry ");

    size_t size = accessions.size();
    size_t last = size - 1;

    for ( size_t i = 0; i < size; ) {
        text << accessions[i];
        ++i;
        if ( i < size ) {
            text << ((i == last) ? " and " : ", ");
        }
    }

    return CNcbiOstrstreamToString(text);
}


string CCommentItem::GetStringForBankIt(const CUser_object& uo)
{
    if ( !uo.CanGetType()  ||  !uo.GetType().IsStr()  ||
         uo.GetType().GetStr() != "Submission" ) {
        return kEmptyStr;
    }

    const string* uvc = 0, *bic = 0;
    if ( uo.HasField("UniVecComment") ) {
        const CUser_field& uf = uo.GetField("UniVecComment");
        if ( uf.CanGetData()  &&  uf.GetData().IsStr() ) {
            uvc = &(uf.GetData().GetStr());
        } 
    }
    if ( uo.HasField("AdditionalComment") ) {
        const CUser_field& uf = uo.GetField("AdditionalComment");
        if ( uf.CanGetData()  &&  uf.GetData().IsStr() ) {
            bic = &(uf.GetData().GetStr());
        } 
    }

    CNcbiOstrstream text;
    if ( uvc != 0  &&  bic != 0 ) {
        text << "Vector Explanation: " << *uvc << "~Bankit Comment: " << *bic;
    } else if ( uvc != 0 ) {
        text << "Vector Explanation: " << *uvc;
    } else if ( bic != 0 ) {
         text << "Bankit Comment: " << *bic;
    }

    return CNcbiOstrstreamToString(text);
}



static void s_GetAssemblyInfo(const CUser_object& uo,
                              string& s,
                              CCommentItem::ECommentFormat format)
{
    s.clear();

    typedef vector < pair<const string*, bool> > TAssemblyInfo;
    TAssemblyInfo assembly;
    if ( uo.HasField("Assembly") ) {
        const CUser_field& field = uo.GetField("Assembly");
        if ( !field.GetData().IsFields() ) {
            return;
        }
        ITERATE (CUser_field::C_Data::TFields, fit,
                 field.GetData().GetFields()) {
            if ( !(*fit)->GetData().IsFields() ) {
                continue;
            }
            ITERATE (CUser_field::C_Data::TFields, it,
                     (*fit)->GetData().GetFields()) {
                const CUser_field& uf = **it;
                if ( !uf.CanGetLabel()  ||  !uf.GetLabel().IsStr() ) {
                    continue;
                }
                const string& label = uf.GetLabel().GetStr();
                if ( label == "accession"  ||  label == "name" ) {
                    bool is_accn = (label == "accession");
                    if ( uf.GetData().IsStr()  &&
                         !uf.GetData().GetStr().empty() ) {
                        assembly.push_back(make_pair(&uf.GetData().GetStr(), is_accn));
                    }
                }
            }
        }
    }
    if ( assembly.size() > 0 ) {
        CNcbiOstrstream oss;
        oss << " The reference sequence was derived from ";
        size_t assembly_size = assembly.size();
        for ( size_t i = 0; i < assembly_size; ++i ) {
            if ( i > 0  ) {
                oss << ((i < assembly_size - 1) ? ", " : " and ");
            }
            const string& acc = *(assembly[i].first);
            NcbiId(oss, acc, format == CCommentItem::eFormat_Html);
        }
        oss << '.';

        s = (string)(CNcbiOstrstreamToString(oss));
    }
}


CCommentItem::TRefTrackStatus CCommentItem::GetRefTrackStatus
(const CUser_object& uo,
 string* st)
{
    TRefTrackStatus retval = eRefTrackStatus_Unknown;
    if ( st != 0 ) {
        st->erase();
    }
    if ( !uo.HasField("Status") ) {
        return retval;
    }

    const CUser_field& field = uo.GetField("Status");
    if ( field.GetData().IsStr() ) {
        string status = field.GetData().GetStr();
        if (NStr::EqualNocase(status, "Inferred")) { 
            retval = eRefTrackStatus_Inferred;
        } else if (NStr::EqualNocase(status, "Provisional")) {
            retval = eRefTrackStatus_Provisional;
        } else if (NStr::EqualNocase(status, "Predicted")) {
            retval = eRefTrackStatus_Predicted;
        } else if (NStr::EqualNocase(status, "Pipeline")) {
            retval = eRefTrackStatus_Pipeline;
        } else if (NStr::EqualNocase(status, "Validated")) {
            retval = eRefTrackStatus_Validated;
        } else if (NStr::EqualNocase(status, "Reviewed")) {
            retval = eRefTrackStatus_Reviewed;
        } else if (NStr::EqualNocase(status, "Model")) {
            retval = eRefTrackStatus_Model;
        } else if (NStr::EqualNocase(status, "WGS")) {
            retval = eRefTrackStatus_WGS;
        }

        if ( st != 0  &&  retval != eRefTrackStatus_Unknown ) {
            *st = NStr::ToUpper(status);
        }
    }

    return retval;
}


string CCommentItem::GetStringForRefTrack
(const CUser_object& uo,
 const CBioseq_Handle& bsh,
 ECommentFormat format)
{
    if ( !uo.IsSetType()  ||  !uo.GetType().IsStr()  ||
         uo.GetType().GetStr() != "RefGeneTracking") {
        return kEmptyStr;
    }

    TRefTrackStatus status = eRefTrackStatus_Unknown;
    string status_str;
    status = GetRefTrackStatus(uo, &status_str);
    if ( status == eRefTrackStatus_Unknown ) {
        return kEmptyStr;
    }

    string collaborator;
    if ( uo.HasField("Collaborator") ) {
        const CUser_field& colab_field = uo.GetField("Collaborator");
        if ( colab_field.GetData().IsStr() ) {
            collaborator = colab_field.GetData().GetStr();
        }
    }

    string source;
    if ( uo.HasField("GenomicSource") ) {
        const CUser_field& source_field = uo.GetField("GenomicSource");
        if ( source_field.GetData().IsStr() ) {
            source = source_field.GetData().GetStr();
        }
    }

    string identical_to;
    if (uo.HasField("IdenticalTo")) {
        const CUser_field& uf = uo.GetField("IdenticalTo");
        ITERATE (CUser_field::TData::TFields, it, uf.GetData().GetFields()) {
            if ( !(*it)->GetData().IsFields() ) {
                continue;
            }
            ITERATE (CUser_field::TData::TFields, i, (**it).GetData().GetFields()) {
                const CUser_field& sub = **i;
                if (sub.GetLabel().GetStr() == "accession") {
                    identical_to = sub.GetData().GetStr();
                    break;
                }
                if (sub.GetLabel().GetStr() == "gi") {
                    identical_to = "gi:" +
                        NStr::IntToString(sub.GetData().GetInt());
                }
            }
        }
    }

    const string *refseq = (format == eFormat_Html ? &kRefSeqLink : &kRefSeq);

    string build_num = CGenomeAnnotComment::GetGenomeBuildNumber(bsh);

    CNcbiOstrstream oss;
    if (status == eRefTrackStatus_Pipeline) {
        oss << *refseq << " INFORMATION: ";
    } else {
        oss << status_str << ' ' << *refseq << ": ";
    }
    switch ( status ) {
    case eRefTrackStatus_Inferred:
        oss << "This record is predicted by genome sequence analysis and is "
            << "not yet supported by experimental evidence.";
        break;
    case eRefTrackStatus_Pipeline:
        if ( !build_num.empty() ) {
            oss << "Features on this sequence have been produced for build "
                << build_num << " of the NCBI's genome annotation"
                << " [see " << "documentation" << "].";
        } else {
            oss << "NCBI contigs are derived from assembled genomic sequence data."
                << "~Also see:~"
                << "    Documentation of NCBI's Annotation Process~ ";
        }
        break;
    case eRefTrackStatus_Provisional:
        if (collaborator.empty()) {
            oss << "This record has not yet been subject to final NCBI review.";
        } else {
            oss << "This record is based on preliminary "
                "annotation provided by " << collaborator << '.';
        }
        break;
    case eRefTrackStatus_Predicted:
        oss << "This record has not been reviewed and the function is unknown.";
        break;
    case eRefTrackStatus_Validated:
        oss << "This record has undergone validation or preliminary review.";
        break;
    case eRefTrackStatus_Reviewed:
        oss << "This record has been curated by " 
            << (collaborator.empty() ? "NCBI staff" : collaborator) << '.';
        break;
    case eRefTrackStatus_Model:
        oss << "This record is predicted by automated computational analysis.";
        break;
    case eRefTrackStatus_WGS:
        oss << "This record is provided to represent a collection of "
            << "whole genome shotgun sequences.";
        break;
    default:
        break;
    }

    if ( status != eRefTrackStatus_Reviewed  &&
         status != eRefTrackStatus_Provisional  &&
         !collaborator.empty() ) {
        oss << " This record has been curated by " << collaborator << '.';
    }

    if ( !source.empty() ) {
        oss << " This record is derived from an annotated genomic sequence ("
            << source << ").";
    }

    if ( !identical_to.empty() ) {
        oss << " The reference sequence is identical to "
            << identical_to << ".";
    }

    {{
         /// add our assembly info
         string s;
         s_GetAssemblyInfo(uo, s, format);
         oss << s;
     }}

    /// check for a concomitant RefSeqGene item
    for (CSeqdesc_CI desc_it(bsh, CSeqdesc::e_User);
         desc_it;  ++desc_it) {
        const CUser_object& obj = desc_it->GetUser();
        if (obj.IsSetType()  &&  obj.GetType().IsStr()  &&
            obj.GetType().GetStr() == "RefSeqGene") {
            CConstRef<CUser_field> f = obj.GetFieldRef("Status");
            if (f  &&  f->GetData().IsStr()) {
                const string& status = f->GetData().GetStr();
                if (status == "Reference Standard") {
                    oss << "~This sequence is a reference standard in "
                        "the RefSeqGene project.";
                }
            }
        }
    }

    return CNcbiOstrstreamToString(oss);
}


string CCommentItem::GetStringForWGS(CBioseqContext& ctx)
{
    static const string default_str = "?";

    if (!ctx.IsWGSMaster()) {
        return kEmptyStr;
    }

    const string& wgsaccn = ctx.GetWGSMasterAccn();
    const string& wgsname = ctx.GetWGSMasterName();

    if (NStr::IsBlank(wgsaccn)  ||  NStr::IsBlank(wgsname)) {
        return kEmptyStr;
    }

    const string* taxname = &default_str;
    for (CSeqdesc_CI it(ctx.GetHandle(), CSeqdesc::e_Source); it; ++it) {
        const CBioSource& src = it->GetSource();
        if (src.IsSetOrg()  &&  src.GetOrg().IsSetTaxname()  &&
            !NStr::IsBlank(src.GetOrg().GetTaxname()) ) {
            taxname = &(src.GetOrg().GetTaxname());
        }
    }

    const string* first = &default_str, *last = &default_str;
    for (CSeqdesc_CI it(ctx.GetHandle(), CSeqdesc::e_User); it; ++it) {
        const CUser_object& uo = it->GetUser();
        if (uo.IsSetType()  &&  uo.GetType().IsStr()  &&
            NStr::EqualNocase(uo.GetType().GetStr(), "WGSProjects")) {
            if (uo.HasField("WGS_accession_first")) {
                const CUser_field& uf = uo.GetField("WGS_accession_first");
                if (uf.IsSetData()  &&  uf.GetData().IsStr()  &&
                    !NStr::IsBlank(uf.GetData().GetStr()) ) {
                    first = &(uf.GetData().GetStr());
                }
            }
            if (uo.HasField("WGS_accession_last")) {
                const CUser_field& uf = uo.GetField("WGS_accession_last");
                if (uf.IsSetData()  &&  uf.GetData().IsStr()  &&
                    !NStr::IsBlank(uf.GetData().GetStr())) {
                    last = &(uf.GetData().GetStr());
                }
            }
        }
    }

    string version = (wgsname.length() == 15) ? 
        wgsname.substr(7, 2) : wgsname.substr(4, 2);

    CNcbiOstrstream text;
    text << "The " << *taxname 
         << " whole genome shotgun (WGS) project has the project accession " 
         << wgsaccn << ".  This version of the project (" << version 
         << ") has the accession number " << wgsname << ",";
    if (*first != *last) {
        text << " and consists of sequences " << *first << "-" << *last;
    } else {
        text << " and consists of sequence " << *first;
    }

    return CNcbiOstrstreamToString(text);
}


string CCommentItem::GetStringForMolinfo(const CMolInfo& mi, CBioseqContext& ctx)
{
    _ASSERT(mi.CanGetCompleteness());

    bool is_prot = ctx.IsProt();

    switch ( mi.GetCompleteness() ) {
    case CMolInfo::eCompleteness_complete:
        return "COMPLETENESS: full length";

    case CMolInfo::eCompleteness_partial:
        return "COMPLETENESS: not full length";

    case CMolInfo::eCompleteness_no_left:
        return (is_prot ? "COMPLETENESS: incomplete on the amino end" :
                          "COMPLETENESS: incomplete on the 5' end");

    case CMolInfo::eCompleteness_no_right:
        return (is_prot ? "COMPLETENESS: incomplete on the carboxy end" :
                          "COMPLETENESS: incomplete on the 3' end");

    case CMolInfo::eCompleteness_no_ends:
        return "COMPLETENESS: incomplete on both ends";

    case CMolInfo::eCompleteness_has_left:
        return (is_prot ? "COMPLETENESS: complete on the amino end" :
                          "COMPLETENESS: complete on the 5' end");

    case CMolInfo::eCompleteness_has_right:
        return (is_prot ? "COMPLETENESS: complete on the carboxy end" :
                          "COMPLETENESS: complete on the 3' end");

    default:
        return "COMPLETENESS: unknown";
    }

    return kEmptyStr;
}


string CCommentItem::GetStringForHTGS(CBioseqContext& ctx)
{
    SDeltaSeqSummary summary;
    if (ctx.IsDelta()) {
        GetDeltaSeqSummary(ctx.GetHandle(), summary);
    }

    CMolInfo::TTech tech = ctx.GetTech();

    CNcbiOstrstream text;

    if ( tech == CMolInfo::eTech_htgs_0 ) {
        if ( summary.num_segs > 0 ) {
            text << "* NOTE: This record contains " << (summary.num_gaps + 1) << " individual~"
                 << "* sequencing reads that have not been assembled into~"
                 << "* contigs. Runs of N are used to separate the reads~"
                 << "* and the order in which they appear is completely~"
                 << "* arbitrary. Low-pass sequence sampling is useful for~"
                 << "* identifying clones that may be gene-rich and allows~"
                 << "* overlap relationships among clones to be deduced.~"
                 << "* However, it should not be assumed that this clone~"
                 << "* will be sequenced to completion. In the event that~"
                 << "* the record is updated, the accession number will~"
                 << "* be preserved.";
        }
        text << "~";
        text << summary.text;
    } else if ( tech == CMolInfo::eTech_htgs_1 ) {
        text << "* NOTE: This is a \"working draft\" sequence.";
        if ( summary.num_segs > 0 ) {
            text << " It currently~"
                 << "* consists of " << (summary.num_gaps + 1) << " contigs. The true order of the pieces~"
                 << "* is not known and their order in this sequence record is~"
                 << "* arbitrary. Gaps between the contigs are represented as~"
                 << "* runs of N, but the exact sizes of the gaps are unknown.";
        }
        text << "~* This record will be updated with the finished sequence~"
             << "* as soon as it is available and the accession number will~"
             << "* be preserved."
             << "~"
             << summary.text;
    } else if ( tech == CMolInfo::eTech_htgs_2 ) {
        text << "* NOTE: This is a \"working draft\" sequence.";
        if ( summary.num_segs > 0 ) {
            text << " It currently~* consists of " << (summary.num_gaps + 1) 
                 << " contigs. Gaps between the contigs~"
                 << "* are represented as runs of N. The order of the pieces~"
                 << "* is believed to be correct as given, however the sizes~"
                 << "* of the gaps between them are based on estimates that have~"
                 << "* provided by the submittor.";
        }
        text << "~* This sequence will be replaced~"
             << "* by the finished sequence as soon as it is available and~"
             << "* the accession number will be preserved."
             << "~"
             << summary.text;
    } else if ( !GetTechString(tech).empty() ) {
        text << "Method: " << GetTechString(tech) << ".";
    }

    string comment = CNcbiOstrstreamToString(text);
    ConvertQuotes(comment);
    ncbi::objects::AddPeriod(comment);

    return comment;
}


string CCommentItem::GetStringForModelEvidance
(const SModelEvidance& me,
 ECommentFormat format)
{
    const string *refseq = (format == eFormat_Html ? &kRefSeqLink : &kRefSeq);

    CNcbiOstrstream text;

    text << "MODEL " << *refseq << ":  " << "This record is predicted by "
         << "automated computational analysis. This record is derived from "
         << "a genomic sequence (" << me.name << ")";
    if ( !me.method.empty() ) {
        text << " annotated using gene prediction method: " << me.method;
    }

    if ( me.mrnaEv  ||  me.estEv ) {
        text << ", supported by ";
        if ( me.mrnaEv  &&  me.estEv ) {
            text << "mRNA and EST ";
        } else if ( me.mrnaEv ) {
            text << "mRNA ";
        } else {
            text << "EST ";
        }
        // !!! for html we need much more !!!
        text << "evidence";
    }

    text << ".~Also see:~"
        << "    Documentation of NCBI's Annotation Process~    ";

    return CNcbiOstrstreamToString(text);
}

static bool s_GetEncodeValues
(string& chromosome,
 string& assembly_date,
 string& ncbi_annotation,
 CBioseqContext& ctx)
{
    _ASSERT(ctx.IsEncode());

    const CUser_object& uo = ctx.GetEncode();
    if (uo.HasField("AssemblyDate")) {
        const CUser_field& ad = uo.GetField("AssemblyDate");
        if (ad.IsSetData()  &&  ad.GetData().IsStr()) {
            assembly_date = ad.GetData().GetStr();
        }
    } else {
        return false;
    }
    if (uo.HasField("NcbiAnnotation")) {
        const CUser_field& na = uo.GetField("NcbiAnnotation");
        if (na.IsSetData()  &&  na.GetData().IsStr()) {
            assembly_date = na.GetData().GetStr();
        }
    } else {
        return false;
    }

    const string* name = NULL;
    for (CSeqdesc_CI it(ctx.GetHandle(), CSeqdesc::e_Source); it; ++it) {
        const CBioSource& bio = it->GetSource();
        ITERATE (CBioSource::TSubtype, st, bio.GetSubtype()) {
            if ((*st)->GetSubtype() == CSubSource::eSubtype_chromosome) {
                name = &(*st)->GetName();
                break;
            }
        }
    }
    if (name != NULL) {
        chromosome = *name;
    } else {
        return false;
    }

    if (NStr::IsBlank(chromosome)) {
        chromosome = "?";
    }
    if (NStr::IsBlank(assembly_date)) {
        assembly_date = "?";
    }
    if (NStr::IsBlank(ncbi_annotation)) {
        ncbi_annotation = "?";
    }
    return true;
}


string CCommentItem::GetStringForEncode(CBioseqContext& ctx)
{
    if (!ctx.IsEncode()) {
        return kEmptyStr;
    }

    CNcbiOstrstream str;
    str << "REFSEQ:  This record was provided by the ENCODE project.";

    string chromosome, assembly_date, ncbi_annotation;
    if (s_GetEncodeValues(chromosome, assembly_date, ncbi_annotation, ctx)) {
        str << "  It is defined by coordinates on the sequence of chromosome"
            << chromosome << "from the" << assembly_date
            << "assembly of the human genome (NCBI build" << ncbi_annotation
            << ").";
    }
    return CNcbiOstrstreamToString(str);
}

/***************************************************************************/
/*                                 PROTECTED                               */
/***************************************************************************/


void CCommentItem::x_GatherInfo(CBioseqContext& ctx)
{
    const CObject* obj = GetObject();
    if (obj == NULL) {
        return;
    }
    const CSeqdesc* desc = dynamic_cast<const CSeqdesc*>(obj);
    if (desc != NULL) {
        x_GatherDescInfo(*desc);
    } else {
        const CSeq_feat* feat = dynamic_cast<const CSeq_feat*>(obj);
        if (feat != NULL) {
            x_GatherFeatInfo(*feat, ctx);
        }
    }
}

// turns data into comment lines (not line-wrapped)
// result in out_lines
// out_prefix_len holds the length of the part up to the space after the double-colon
static 
void s_GetStrForStructuredComment( 
    const CUser_object::TData &data, 
    list<string> &out_lines,
    int &out_prefix_len,
    const bool is_first )
{
    // default prefix and suffix
    const char* prefix = "##Metadata-START##";
    const char* suffix = "##Metadata-END##";

    // First, figure out the longest label so we know how to format it
    // (and set the prefix and suffix while we're at it)
    string::size_type longest_label_len = 1;
    ITERATE( CUser_object::TData, it_for_len, data ) {
        if( (*it_for_len)->GetLabel().IsStr() && 
                (*it_for_len)->GetData().IsStr() && ! (*it_for_len)->GetData().GetStr().empty() ) {
            const string &label = (*it_for_len)->GetLabel().GetStr();

            if( label == "StructuredCommentPrefix" ) {
                prefix = (*it_for_len)->GetData().GetStr().c_str();
            } else if( label == "StructuredCommentSuffix" ) {
                suffix = (*it_for_len)->GetData().GetStr().c_str();
            } else {
                const string::size_type label_len = label.length();
                if( (label_len > longest_label_len) && (label_len <= 45) ) {
                    longest_label_len = label_len;
                }
            }
        }
    }
    out_prefix_len = (longest_label_len + 4); // "+4" because we add " :: " after the prefix

    if( ! is_first ) {
        out_lines.push_back( "\n" );
    }
    out_lines.push_back( prefix );
    out_lines.back().append( "\n" );

    ITERATE( CUser_object::TData, it, data ) {
        
        // skip if no label
        if( ! (*it)->GetLabel().IsStr() || (*it)->GetLabel().GetStr().empty() ) {
            continue;
        }

        // skip if no data
        if( ! (*it)->GetData().IsStr() || (*it)->GetData().GetStr().empty() ) {
            continue;
        }

        // special fields are skipped
        if( (*it)->GetLabel().GetStr() == "StructuredCommentPrefix" || 
                (*it)->GetLabel().GetStr() == "StructuredCommentSuffix" ) {
            continue;
        }

        // crate the next line that we're going to set the contents of
        out_lines.push_back( (*it)->GetLabel().GetStr() );
        string &next_line = out_lines.back();

        next_line.resize( max( next_line.size(), longest_label_len), ' ' );
        next_line.append( " :: " );
        next_line.append( (*it)->GetData().GetStr() );
        next_line.append( "\n" );

        ExpandTildes(next_line, eTilde_comment);
    }

    out_lines.push_back( suffix );
    out_lines.back().append( "\n" );
}

void CCommentItem::x_GatherDescInfo(const CSeqdesc& desc)
{
    // true for most desc infos
    EPeriod can_add_period = ePeriod_Add;

    string prefix, str, suffix;
    switch ( desc.Which() ) {
    case CSeqdesc::e_Comment:
        {{
            if (!NStr::IsBlank(desc.GetComment())) {
                str = desc.GetComment();
                TrimSpacesAndJunkFromEnds(str);
                ConvertQuotes(str);
                ncbi::objects::AddPeriod(str);
            }
        }}
        break;

    case CSeqdesc::e_Maploc:
        {{
            const CDbtag& dbtag = desc.GetMaploc();
            if ( dbtag.CanGetTag() ) {
                const CObject_id& oid = dbtag.GetTag();
                if ( oid.IsStr() ) {
                    prefix = "Map location: ";
                    str = oid.GetStr();
                } else if ( oid.IsId()  &&  dbtag.CanGetDb() ) {
                    prefix = "Map location: (Database ";
                    str = dbtag.GetDb();
                    suffix = "; id # " + NStr::IntToString(oid.GetId()) + ").";
                }
            }
        }}
        break;

    case CSeqdesc::e_Region:
        {{
            prefix = "Region: ";
            str = desc.GetRegion();
            NStr::ReplaceInPlace(str, "\"", "\'");
            ncbi::objects::AddPeriod(str);
        }}
        break;

    case CSeqdesc::e_Name:
        {{
            prefix = "Name: ";
            str = desc.GetName();
            ncbi::objects::AddPeriod(str);
        }}
        break;

    case CSeqdesc::e_User:
        {{
            const CSeqdesc_Base::TUser &userObject = desc.GetUser();

            // make sure the user object is really of type StructuredComment
            const CUser_object::TType &type = userObject.GetType();
            if( type.IsStr() && type.GetStr() == "StructuredComment" ) {
                s_GetStrForStructuredComment( userObject.GetData(),  
                    m_Comment, m_CommentInternalIndent, IsFirst() );
                SetNeedPeriod( false );
                can_add_period = ePeriod_NoAdd;
                return; // special case because multiple lines
            }
        }}
        break;

    default:
        break;
    }

    if (str.empty()  ||  str == ".") {
        return;
    }
    x_SetCommentWithURLlinks(prefix, str, suffix, can_add_period);
    
}


void CCommentItem::x_GatherFeatInfo(const CSeq_feat& feat, CBioseqContext& ctx)
{
    if (!feat.GetData().IsComment()  ||
        !feat.CanGetComment()        ||
        NStr::IsBlank(feat.GetComment())) {
        return;
    }

    x_SetCommentWithURLlinks(kEmptyStr, feat.GetComment(), kEmptyStr, ePeriod_Add);
}


void CCommentItem::x_SetSkip(void)
{
    CFlatItem::x_SetSkip();
    swap(m_First, sm_FirstComment);
}


void CCommentItem::x_SetComment(const string& comment)
{
    m_Comment.clear();
    m_Comment.push_back( comment );
    ExpandTildes(m_Comment.back(), eTilde_comment);;
}


void CCommentItem::x_SetCommentWithURLlinks
(const string& prefix,
 const string& str,
 const string& suffix,
 EPeriod can_add_period)
{
    // !!! test for html - find links within the comment string
    string comment = prefix;
    comment += str;
    comment += suffix;

    ExpandTildes(comment, eTilde_comment);
    if (NStr::IsBlank(comment)) {
        return;
    }

    if( can_add_period == ePeriod_Add ) {
        size_t pos = comment.find_last_not_of(" \n\t\r.~");
        if (pos != comment.length() - 1) {
            size_t period = comment.find_last_of('.');
            bool add_period = period > pos;
            if (add_period  &&  !NStr::EndsWith(str, "...")) {
                ncbi::objects::AddPeriod(comment);
            }
        }
    }
    
    ConvertQuotes( comment );

    m_Comment.clear();
    m_Comment.push_back( comment );
}

bool CCommentItem::x_IsCommentEmpty(void) const
{
    ITERATE(list<string>, it, m_Comment) {
        if( ! m_Comment.empty() ) {
            return false;
        }
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////
//
// Derived Classes

// --- CGenomeAnnotComment

CGenomeAnnotComment::CGenomeAnnotComment
(CBioseqContext& ctx,
 const string& build_num) :
    CCommentItem(ctx), m_GenomeBuildNumber(build_num)
{
    x_GatherInfo(ctx);
}


string CGenomeAnnotComment::GetGenomeBuildNumber(const CUser_object& uo)
{
    if ( uo.IsSetType()  &&  uo.GetType().IsStr()  &&
         uo.GetType().GetStr() == "GenomeBuild" ) {
        if ( uo.HasField("NcbiAnnotation") ) {
            string build_num;
            const CUser_field& uf = uo.GetField("NcbiAnnotation");
            if ( uf.CanGetData()  &&  uf.GetData().IsStr()  &&
                 !uf.GetData().GetStr().empty() ) {
                build_num = uf.GetData().GetStr();
            }

            if ( uo.HasField("NcbiVersion") ) {
                const CUser_field& uf = uo.GetField("NcbiVersion");
                if ( uf.CanGetData()  &&  uf.GetData().IsStr()  &&
                     !uf.GetData().GetStr().empty() ) {
                    build_num += " version ";
                    build_num += uf.GetData().GetStr();
                }
            }
            return build_num;

        } else if ( uo.HasField("Annotation") ) {
            const CUser_field& uf = uo.GetField("Annotation");
            if ( uf.CanGetData()  &&  uf.GetData().IsStr()  &&
                 !uf.GetData().GetStr().empty() ) {
                static const string prefix = "NCBI build ";
                if ( NStr::StartsWith(uf.GetData().GetStr(), prefix) ) {
                    return uf.GetData().GetStr().substr(prefix.length());
                }
            }
        }
    }
    return kEmptyStr;
}


string CGenomeAnnotComment::GetGenomeBuildNumber(const CBioseq_Handle& bsh)
{
    for (CSeqdesc_CI it(bsh, CSeqdesc::e_User);  it;  ++it) {
        const CUser_object& uo = it->GetUser();
        string s = GetGenomeBuildNumber(uo);
        if ( !s.empty() ) {
            return s;
        }
    }

    return kEmptyStr;
}

void CGenomeAnnotComment::x_GatherInfo(CBioseqContext& ctx)
{
    const string *refseq = ctx.Config().DoHTML() ? &kRefSeqLink : &kRefSeq;

    CNcbiOstrstream text;

    text << "GENOME ANNOTATION " << *refseq << ": ";
    if ( !m_GenomeBuildNumber.empty() ) {
         text << "Features on this sequence have been produced for build "
              << m_GenomeBuildNumber << " of the NCBI's genome annotation"
              << " [see " << "documentation" << "].";
    } else {
        text << "NCBI contigs are derived from assembled genomic sequence data."
             << "~Also see:~"
             << "    Documentation of NCBI's Annotation Process~ ";
    }

    /// add our assembly info
    for (CSeqdesc_CI desc_it(ctx.GetHandle(), CSeqdesc::e_User);
         desc_it;  ++desc_it) {
        const CUser_object& uo = desc_it->GetUser();
        if ( !uo.IsSetType()  ||  !uo.GetType().IsStr()  ||
             uo.GetType().GetStr() != "RefGeneTracking") {
            continue;
        }

        string s;
        s_GetAssemblyInfo(uo, s,
                          ctx.Config().DoHTML() ?
                          CCommentItem::eFormat_Html :
                          CCommentItem::eFormat_Text);
        text << s;
        break;
    }

    string s = (string)(CNcbiOstrstreamToString(text));
    x_SetComment(s);
}


// --- CHistComment

CHistComment::CHistComment
(EType type,
 const CSeq_hist& hist,
 CBioseqContext& ctx) : 
    CCommentItem(ctx), m_Type(type), m_Hist(&hist)
{
    x_GatherInfo(ctx);
    m_Hist.Reset();
}


string s_CreateHistCommentString
(const string& prefix,
 const string& suffix,
 const CSeq_hist_rec& hist,
 bool do_html)
{
    //if (!hist.CanGetDate()  ||  !hist.CanGetIds()) {
    //    return "???";
    //}

    string date;
    if (hist.IsSetDate()) {
        hist.GetDate().GetDate(&date, "%{%3N%|???%} %{%D%|??%}, %{%4Y%|????%}");
    }

    vector<int> gis;
    ITERATE (CSeq_hist_rec::TIds, id, hist.GetIds()) {
        if ( (*id)->IsGi() ) {
            gis.push_back((*id)->GetGi());
        }
    }

    CNcbiOstrstream text;

    text << prefix << ((gis.size() > 1) ? " or before " : " ") << date 
         << ' ' << suffix;

    if ( gis.empty() ) {
        text << " gi:?";
        return CNcbiOstrstreamToString(text);
    }

    for ( size_t count = 0; count < gis.size(); ++count ) {
        if ( count != 0 ) {
            text << ",";
        }
        text << " gi:";
        NcbiId(text, gis[count], do_html);
    }
    text << '.' << endl;

    return CNcbiOstrstreamToString(text);
}

void CHistComment::x_GatherInfo(CBioseqContext& ctx)
{
    _ASSERT(m_Hist);

    switch ( m_Type ) {
    case eReplaced_by:
        if( ctx.IsWGSMaster() ) {
            x_SetComment(s_CreateHistCommentString(
                "[WARNING] On",
                "this project was updated. The new version is",
                m_Hist->GetReplaced_by(),
                ctx.Config().DoHTML()));
        } else {
            x_SetComment(s_CreateHistCommentString(
                "[WARNING] On",
                "this sequence was replaced by",
                m_Hist->GetReplaced_by(),
                ctx.Config().DoHTML()));
        }
        break;
    case eReplaces:
        x_SetComment(s_CreateHistCommentString(
            "On",
            "this sequence version replaced",
            m_Hist->GetReplaces(),
            ctx.Config().DoHTML()));
        break;
    }
}


// --- CGsdbComment

CGsdbComment::CGsdbComment(const CDbtag& dbtag, CBioseqContext& ctx) :
    CCommentItem(ctx), m_Dbtag(&dbtag)
{
    x_GatherInfo(ctx);
}


void CGsdbComment::x_GatherInfo(CBioseqContext&)
{
    if (m_Dbtag->IsSetTag()  &&  m_Dbtag->GetTag().IsId()) {
        string id = NStr::IntToString(m_Dbtag->GetTag().GetId());
        x_SetComment("GSDB:S:" + id);
    } else {
        x_SetSkip();
    }
}


// --- CLocalIdComment

CLocalIdComment::CLocalIdComment(const CObject_id& oid, CBioseqContext& ctx) :
    CCommentItem(ctx, false), m_Oid(&oid)
{
    x_GatherInfo(ctx);
}


void CLocalIdComment::x_GatherInfo(CBioseqContext&)
{
    CNcbiOstrstream msg;

    switch ( m_Oid->Which() ) {
    case CObject_id::e_Id:
        msg << "LocalID: " << m_Oid->GetId();    
        break;
    case CObject_id::e_Str:
        if ( m_Oid->GetStr().length() < 1000 ) {
            msg << "LocalID: " << m_Oid->GetStr();
        } else {
            msg << "LocalID string too large";
        }
        break;
    default:
        break;
    }
    x_SetComment(CNcbiOstrstreamToString(msg));
}

END_SCOPE(objects)
END_NCBI_SCOPE
