/* $Id: cleanup_gbqual.cpp 240004 2011-02-02 16:24:02Z rafanovi $
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
 * Author:  Robert G. Smith
 *
 * File Description:
 *   Implementation of BasicCleanup for GBQuals in Seq-feat's.
 *
 */

#include <ncbi_pch.hpp>
#include "cleanup_utils.hpp"
#include <objects/general/Object_id.hpp>
#include <objects/general/Dbtag.hpp>
#include <objects/seqfeat/Cdregion.hpp>
#include <objects/seqfeat/Genetic_code.hpp>
#include <objects/seqfeat/Genetic_code_table.hpp>
#include <objects/seqfeat/Gb_qual.hpp>
#include <objects/seqfeat/Imp_feat.hpp>
#include <objects/seqfeat/Seq_feat.hpp>
#include <objects/seqfeat/SeqFeatXref.hpp>
#include <objects/seqfeat/Trna_ext.hpp>
#include <objects/seqfeat/Code_break.hpp>
#include <util/static_map.hpp>

#include <objects/seqfeat/RNA_ref.hpp>
#include <objects/seqloc/Seq_loc.hpp>
#include <objects/general/User_object.hpp>
#include <objects/general/User_field.hpp>
#include <objects/seq/seqport_util.hpp>
#include <objects/misc/sequence_macros.hpp>
#include <vector>

#include <objmgr/util/seq_loc_util.hpp>

#include <objmgr/feat_ci.hpp>
#include <objmgr/annot_ci.hpp>

#include "cleanupp.hpp"


BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

// seq-feat.qual


bool CCleanup_imp::BasicCleanup(CGene_ref& gene, const CGb_qual& gb_qual)
{
    const string& qual = gb_qual.GetQual();
    const string& val  = gb_qual.GetVal();

    bool retval = false;
    if (NStr::EqualNocase(qual, "map")) {
        if (! (gene.IsSetMaploc()  ||  NStr::IsBlank(val)) ) {
            retval = true;
            gene.SetMaploc(val);
        }
    } else if (NStr::EqualNocase(qual, "allele")) {
        if ( ! (gene.IsSetAllele()  ||  NStr::IsBlank(val)) ) {
            retval = true;
            gene.SetAllele(val);
        }
    } else if (NStr::EqualNocase(qual, "locus_tag")) {
        if ( ! (gene.IsSetLocus_tag()  ||  NStr::IsBlank(val)) ) {
            retval = true;
            gene.SetLocus_tag(val);
        }
    } else if (NStr::EqualNocase(qual, "gene")  &&  ! NStr::IsBlank(val)) {
        retval = true;
        if ( ! gene.IsSetLocus() ) {
            gene.SetLocus(val);
        } else if (gene.GetLocus() != val) {
            CGene_ref::TSyn::const_iterator syn_it = 
                find(gene.GetSyn().begin(), gene.GetSyn().end(), val);
            if (syn_it == gene.GetSyn().end()) {
                gene.SetSyn().push_back(val);
            }            
        }
    }
    if (retval) {
        ChangeMade(CCleanupChange::eChangeQualifiers);
    }

    return retval;
}


bool CCleanup_imp::x_ParseCodeBreak(const CSeq_feat& feat, CCdregion& cds, const string& str)
{
    string::size_type aa_pos = NStr::Find(str, "aa:");
    string::size_type len = 0;
    string::size_type loc_pos, end_pos;
    char protein_letter = 'X';
    CRef<CSeq_loc> break_loc;
    
    if (aa_pos == string::npos) {
        aa_pos = NStr::Find (str, ",");
        if (aa_pos != string::npos) {
            aa_pos = NStr::Find (str, ":", aa_pos);
        }
        if (aa_pos != string::npos) {
            aa_pos ++;
        }
    } else {
        aa_pos += 3;
    }

    if (aa_pos != string::npos) {    
        while (aa_pos < str.length() && isspace (str[aa_pos])) {
            aa_pos++;
        }
        while (aa_pos + len < str.length() && isalpha (str[aa_pos + len])) {
            len++;
        }
        if (len != 0) {    
            protein_letter = ValidAminoAcid(str.substr(aa_pos, len));
        }
    }
    
    loc_pos = NStr::Find (str, "(pos:");
    if (loc_pos == string::npos) {
        return false;
    }
    loc_pos += 5;
    while (loc_pos < str.length() && isspace (str[loc_pos])) {
        loc_pos++;
    }
    end_pos = NStr::Find (str, ",", loc_pos);
    if (end_pos == string::npos) {
        break_loc = ReadLocFromText (str.substr(loc_pos), feat.GetLocation().GetId(), m_Scope);
    } else {
        break_loc = ReadLocFromText (str.substr(loc_pos, end_pos - loc_pos), feat.GetLocation().GetId(), m_Scope);
    }
    
    if (break_loc == NULL 
        || sequence::Compare (*break_loc, feat.GetLocation(), m_Scope) != sequence::eContained
        || (break_loc->IsInt() && sequence::GetLength(*break_loc, m_Scope) != 3)) {
        return false;
    }
    
    // need to build code break object and add it to coding region
    CRef<CCode_break> newCodeBreak(new CCode_break());
    CCode_break::TAa& aa = newCodeBreak->SetAa();
    aa.SetNcbieaa(protein_letter);
    newCodeBreak->SetLoc (*break_loc);

    CCdregion::TCode_break& orig_list = cds.SetCode_break();
    orig_list.push_back(newCodeBreak);
    
    ChangeMade(CCleanupChange::eChangeCodeBreak);
    
    return true;
    
}


bool CCleanup_imp::BasicCleanup(CSeq_feat& feat, CCdregion& cds, const CGb_qual& gb_qual)
{
    const string& qual = gb_qual.GetQual();
    const string& val  = gb_qual.GetVal();
    
    // transl_except qual -> Cdregion.code_break
    if (NStr::EqualNocase(qual, "transl_except")) {
        
        return x_ParseCodeBreak(feat, cds, val);
    }

    // codon_start qual -> Cdregion.frame
    if (NStr::EqualNocase(qual, "codon_start")) {
        CCdregion::TFrame frame = cds.GetFrame();
        CCdregion::TFrame new_frame = CCdregion::TFrame(NStr::StringToNumeric(val));
        if (new_frame == CCdregion::eFrame_one  ||
            new_frame == CCdregion::eFrame_two  ||
            new_frame == CCdregion::eFrame_three) {
            if (frame == CCdregion::eFrame_not_set  ||
                (feat.IsSetPseudo()  &&  feat.GetPseudo()  &&  !feat.IsSetProduct())) {
                cds.SetFrame(new_frame);
                ChangeMade(CCleanupChange::eChangeQualifiers);
            }
            return true;
        }
    }

    // transl_table qual -> Cdregion.code
    if (NStr::EqualNocase(qual, "transl_table")) {
        if (cds.IsSetCode()) {
            const CCdregion::TCode& code = cds.GetCode();
            int transl_table = 1;
            ITERATE (CCdregion::TCode::Tdata, it, code.Get()) {
                if ((*it)->IsId()  &&  (*it)->GetId() != 0) {
                    transl_table = (*it)->GetId();
                    break;
                }
            }
            
            if (NStr::EqualNocase(NStr::UIntToString(transl_table), val)) {
                return true;
            }
        } else {
            int new_val = NStr::StringToNumeric(val);
            if (new_val > 0) {
                CRef<CGenetic_code::C_E> gc(new CGenetic_code::C_E);
                gc->SetId(new_val);
                cds.SetCode().Set().push_back(gc);
                
                ChangeMade(CCleanupChange::eChangeGeneticCode);
                return true;
            }
        }
    }

    // look for qualifiers that should be applied to protein feature
    // note - this should be moved to the "indexed" portion of basic cleanup,
    // because it needs to locate another sequence and feature
    if (NStr::Equal(qual, "product") || NStr::Equal (qual, "function") || NStr::Equal (qual, "EC_number")
        || NStr::Equal(qual, "activity") || NStr::Equal (qual, "prot_note")) {

        if (feat.IsSetProduct()) {
            // get protein sequence for product
            CBioseq_Handle prot = m_Scope->GetBioseqHandle(feat.GetProduct());
            if (prot) {
                // replacement prot feature
                CRef<CSeq_feat> prot_feat(new CSeq_feat());

                // find main protein feature
                SAnnotSelector sel(CSeqFeatData::eSubtype_prot);

                CFeat_CI feat_ci (prot, sel);

                if (feat_ci) {            
                    prot_feat->Assign(feat_ci->GetOriginalFeature());

                    bool change_made = false;


                    if (NStr::Equal(qual, "prot_note")) {
                        if (!prot_feat->IsSetComment() || NStr::IsBlank (prot_feat->GetComment())) {
                            prot_feat->SetComment (val);
                        } else {
                            prot_feat->SetComment (prot_feat->GetComment() + "; " + val);
                        }
                        ChangeMade (CCleanupChange::eChangeComment);
                        change_made = true;
                    } else {
                        change_made = BasicCleanup (prot_feat->SetData().SetProt(), gb_qual);
                    }

                    if (change_made) {
                        CSeq_feat_EditHandle efh(feat_ci->GetSeq_feat_Handle());
                        efh.Replace(*prot_feat);
                        return true;
                    }
                }
            }
        } else if (!NStr::Equal(qual, "prot_note")) {
            bool found = false;
            bool change_made = false;
            // find or create prot xref for feature
            EDIT_EACH_SEQFEATXREF_ON_SEQFEAT (it, feat) {
                if ((*it)->IsSetData() && (*it)->GetData().IsProt()) {
                    found = true;
                    change_made = BasicCleanup ((*it)->SetData().SetProt(), gb_qual);
                }
            }
            if (!found) {
                change_made = BasicCleanup (feat.SetProtXref(), gb_qual);
                ChangeMade (CCleanupChange::eAddProtXref);                  
            }
            if (change_made) {
                return true;
            }
        }
    }

    /*if (NStr::EqualNocase(qual, "translation")) {
        return TRUE;
    }*/
    return false;
}


typedef pair <const char *, const int> TTrnaKey;

static const TTrnaKey trna_key_to_subtype [] = {
    TTrnaKey ( "Ala",            'A' ),
    TTrnaKey ( "Alanine",        'A' ),
    TTrnaKey ( "Arg",            'R' ),
    TTrnaKey ( "Arginine",       'R' ),
    TTrnaKey ( "Asn",            'N' ),
    TTrnaKey ( "Asp",            'D' ),
    TTrnaKey ( "Asp or Asn",     'B' ),
    TTrnaKey ( "Asparagine",     'N' ),
    TTrnaKey ( "Aspartate",      'D' ),
    TTrnaKey ( "Aspartic Acid",  'D' ),
    TTrnaKey ( "Asx",            'B' ),
    TTrnaKey ( "Cys",            'C' ),
    TTrnaKey ( "Cysteine",       'C' ),
    TTrnaKey ( "fMet",           'M' ),
    TTrnaKey ( "Gln",            'Q' ),
    TTrnaKey ( "Glu",            'E' ),
    TTrnaKey ( "Glu or Gln",     'Z' ),
    TTrnaKey ( "Glutamate",      'E' ),
    TTrnaKey ( "Glutamic Acid",  'E' ),
    TTrnaKey ( "Glutamine",      'Q' ),
    TTrnaKey ( "Glx",            'Z' ),
    TTrnaKey ( "Gly",            'G' ),
    TTrnaKey ( "Glycine",        'G' ),
    TTrnaKey ( "His",            'H' ),
    TTrnaKey ( "Histidine",      'H' ),
    TTrnaKey ( "Ile",            'I' ),
    TTrnaKey ( "Isoleucine",     'I' ),
    TTrnaKey ( "Leu",            'L' ),
    TTrnaKey ( "Leu or Ile",     'J' ),
    TTrnaKey ( "Leucine",        'L' ),
    TTrnaKey ( "Lys",            'K' ),
    TTrnaKey ( "Lysine",         'K' ),
    TTrnaKey ( "Met",            'M' ),
    TTrnaKey ( "Methionine",     'M' ),
    TTrnaKey ( "OTHER",          'X' ),
    TTrnaKey ( "Phe",            'F' ),
    TTrnaKey ( "Phenylalanine",  'F' ),
    TTrnaKey ( "Pro",            'P' ),
    TTrnaKey ( "Proline",        'P' ),
    TTrnaKey ( "Pyl",            'O' ),
    TTrnaKey ( "Pyrrolysine",    'O' ),
    TTrnaKey ( "Sec",            'U' ),
    TTrnaKey ( "Selenocysteine", 'U' ),
    TTrnaKey ( "Ser",            'S' ),
    TTrnaKey ( "Serine",         'S' ),
    TTrnaKey ( "Ter",            '*' ),
    TTrnaKey ( "TERM",           '*' ),
    TTrnaKey ( "Termination",    '*' ),
    TTrnaKey ( "Thr",            'T' ),
    TTrnaKey ( "Threonine",      'T' ),
    TTrnaKey ( "Trp",            'W' ),
    TTrnaKey ( "Tryptophan",     'W' ),
    TTrnaKey ( "Tyr",            'Y' ),
    TTrnaKey ( "Tyrosine",       'Y' ),
    TTrnaKey ( "Val",            'V' ),
    TTrnaKey ( "Valine",         'V' ),
    TTrnaKey ( "Xle",            'J' ),
    TTrnaKey ( "Xxx",            'X' )
};

typedef CStaticArrayMap <const char*, const int, PNocase_CStr> TTrnaMap;
DEFINE_STATIC_ARRAY_MAP(TTrnaMap, sm_TrnaKeys, trna_key_to_subtype);


static CRef<CTrna_ext> s_ParseTRnaFromAnticodonString (string str, CSeq_feat& feat, CScope *scope)
{
    CRef<CTrna_ext> trna (new CTrna_ext());
    
    if (NStr::IsBlank (str)) return trna;

    if (NStr::StartsWith (str, "(pos:")) {
        string::size_type pos_end = NStr::Find (str, ")");
        if (pos_end != string::npos) {
            string pos_str = str.substr (5, pos_end - 5);
            string::size_type aa_start = NStr::FindNoCase (pos_str, "aa:");
            if (aa_start != string::npos) {
                string abbrev = pos_str.substr (aa_start + 3);
                TTrnaMap::const_iterator t_iter = sm_TrnaKeys.find (abbrev.c_str ());
                if (t_iter == sm_TrnaKeys.end ()) {
                    // unable to parse
                    return trna;
                }
                CRef<CTrna_ext::TAa> aa(new CTrna_ext::TAa);
                aa->SetIupacaa (t_iter->second);
                trna->SetAa(*aa);
                pos_str = pos_str.substr (0, aa_start);
                NStr::TruncateSpacesInPlace (pos_str);
                if (NStr::EndsWith (pos_str, ",")) {
                    pos_str = pos_str.substr (0, pos_str.length() - 1);
                }
            }
            CRef<CSeq_loc> anticodon = ReadLocFromText (pos_str, feat.GetLocation().GetId(), scope);
            if (anticodon == NULL) {
                trna->ResetAa();
            } else {
                trna->SetAnticodon(*anticodon);
            }
        }
    }
    return trna;        
}


static CRef<CTrna_ext> s_ParseTRnaString (string &str)
{
    CRef<CTrna_ext> trna (new CTrna_ext());
    
    if (NStr::IsBlank (str)) return trna;

    /* try the whole string for product */
    TTrnaMap::const_iterator t_iter = sm_TrnaKeys.find (str.c_str ());
    if (t_iter != sm_TrnaKeys.end ()) {
        CRef<CTrna_ext::TAa> aa(new CTrna_ext::TAa);
        aa->SetIupacaa (t_iter->second);
        trna->SetAa(*aa);
        str = "";
        return trna;
    }
    
    if (!NStr::StartsWith (str, "tRNA-")) {
        return trna;
    } else {
        str = str.substr (5);
        string::size_type aa_end = NStr::Find(str, "(");
        

        string abbrev = "";
        if (aa_end == string::npos) {
            abbrev = str;
        } else {
            abbrev = str.substr (0, aa_end);
        }
        NStr::TruncateSpacesInPlace (abbrev);
        TTrnaMap::const_iterator t_iter = sm_TrnaKeys.find (abbrev.c_str ());
        if (t_iter == sm_TrnaKeys.end ()) {
            // couldn't find abbreviation in list
            return trna;
        }
        CRef<CTrna_ext::TAa> aa(new CTrna_ext::TAa);
        aa->SetIupacaa (t_iter->second);

        trna->SetAa( *aa);
        if (aa_end == string::npos) {
            // abbreviation was all there was to find
            str = "";
            return trna;
        }
        // continue parsing
        str = str.substr(aa_end + 1);
        
        string::size_type codons_end = NStr::Find (str, ")");
        if (codons_end != string::npos) {
            string codon_list = str.substr (0, codons_end);
            vector<string> codons;
            vector<int> codon_vals;
            NStr::Tokenize(codon_list, ",", codons);
            bool codons_ok = true;
            for (unsigned int i = 0; i < codons.size() && codons_ok; i++) {
                NStr::TruncateSpacesInPlace (codons[i]);
                if (codons[i].length() != 3) {
                    codons_ok = false;
                } else {
                    NStr::ToUpper (codons[i]);
                    NStr::ReplaceInPlace (codons[i], "T", "U");
                    for (int j = 0; j < 3; j++) {
                        string letter = codons[i].substr(j, 1);
                        if (!NStr::Equal (letter, "A")
                            && !NStr::Equal (letter, "U")
                            && !NStr::Equal (letter, "G")
                            && !NStr::Equal (letter, "C")) {
                            codons_ok = false;
                        }
                    }
                    if (codons_ok) {
                        int residue = CGen_code_table::CodonToIndex (codons[i]);
                        if (residue < 0) {
                            codons_ok = false;
                        } else {
                            codon_vals.push_back (residue);
                        }
                    }
                }
            }
            if (codons_ok) {
                CTrna_ext::TCodon& real_codons = trna->SetCodon ();
                for (unsigned int i = 0; i < codon_vals.size(); i++) {
                    real_codons.push_back (codon_vals[i]);
                }
                str = str.substr (codons_end + 1);
            }
        }
    }
    return trna;
}


bool CCleanup_imp::BasicCleanup(CSeq_feat& feat, CRNA_ref& rna, CGb_qual& gb_qual)
{
    const string& qual = gb_qual.GetQual();

    bool is_std_name = NStr::EqualNocase(qual, "standard_name");
    if (NStr::EqualNocase(qual, "product")  ||  (is_std_name  &&  ! (m_Mode == eCleanup_EMBL  ||  m_Mode == eCleanup_DDBJ) )) {
        if (!gb_qual.IsSetVal()) {
            return false;
        }
        if (rna.IsSetType()) {
            if (rna.GetType() == CRNA_ref::eType_unknown) {
                rna.SetType(CRNA_ref::eType_other);
                ChangeMade(CCleanupChange::eChangeKeywords);
            }
        } else {
            rna.SetType(CRNA_ref::eType_other);
            ChangeMade(CCleanupChange::eChangeKeywords);
        }
        _ASSERT(rna.IsSetType());

        CRNA_ref::TType type = rna.GetType();
        
        if (type == CRNA_ref::eType_other  &&  is_std_name) {
            return false;
        }

        if (type == CRNA_ref::eType_other && rna.IsSetExt() && rna.GetExt().IsName()
            && NStr::Equal (rna.GetExt().GetName(), "misc_RNA")
            && NStr::Equal(qual, "product")) {
            if (NStr::EqualNocase(gb_qual.GetVal(), "its1")
                || NStr::EqualNocase(gb_qual.GetVal(), "its 1")
                || NStr::EqualNocase(gb_qual.GetVal(), "Ribosomal DNA internal transcribed spacer 1")
                || NStr::EqualNocase(gb_qual.GetVal(), "internal transcribed spacer 1 (ITS1)")) {
                gb_qual.SetVal("internal transcribed spacer 1");
                ChangeMade(CCleanupChange::eChangeQualifiers);
                return false;
            } else if (NStr::EqualNocase(gb_qual.GetVal(), "its2")
                || NStr::EqualNocase(gb_qual.GetVal(), "its 2")
                || NStr::EqualNocase(gb_qual.GetVal(), "Ribosomal DNA internal transcribed spacer 2")
                || NStr::EqualNocase(gb_qual.GetVal(), "internal transcribed spacer 2 (ITS2)")) {
                gb_qual.SetVal("internal transcribed spacer 2");
                ChangeMade(CCleanupChange::eChangeQualifiers);
                return false;
            } else if (NStr::EqualNocase(gb_qual.GetVal(), "its3")
                || NStr::EqualNocase(gb_qual.GetVal(), "its 3")
                || NStr::EqualNocase(gb_qual.GetVal(), "Ribosomal DNA internal transcribed spacer 3")
                || NStr::EqualNocase(gb_qual.GetVal(), "internal transcribed spacer 3 (ITS3)")) {
                gb_qual.SetVal("internal transcribed spacer 3");
                ChangeMade(CCleanupChange::eChangeQualifiers);
                return false;
            }
        }


        if (type == CRNA_ref::eType_tRNA) {
            if (rna.IsSetExt()) {
                if (rna.GetExt().IsName()) { 
                    string comment = rna.GetExt().GetName();
                    CRef<CTrna_ext> trna_from_name = s_ParseTRnaString (comment);
                    if (trna_from_name->IsSetAa() && NStr::IsBlank (comment)) {
                        rna.SetExt().SetTRNA (*trna_from_name);
                        ChangeMade (CCleanupChange::eChange_tRna);
                        if (trna_from_name->GetAa().GetIupacaa() == 77 
                            && NStr::Find (rna.GetExt().GetName(), "fMet") != string::npos) {
                            if (!feat.IsSetComment()) {
                                feat.SetComment ("fMet");
                            } else if (NStr::Find (feat.GetComment(), "fMet") == string::npos) {
                                feat.SetComment (feat.GetComment() + "; fMet");
                            }
                            ChangeMade (CCleanupChange::eChangeComment);
                        }
                    }
                }
            } else {                
                string comment = gb_qual.GetVal();
                CRef<CTrna_ext> trna_from_name = s_ParseTRnaString (comment);
                if (trna_from_name->IsSetAa() && NStr::IsBlank (comment)) {
                    rna.SetExt().SetTRNA (*trna_from_name);
                    ChangeMade (CCleanupChange::eChange_tRna);
                    if (trna_from_name->GetAa().GetIupacaa() == 77 
                        && NStr::Find (gb_qual.GetVal(), "fMet") != string::npos) {
                        if (!feat.IsSetComment()) {
                            feat.SetComment ("fMet");
                        } else if (NStr::Find (feat.GetComment(), "fMet") == string::npos) {
                            feat.SetComment (feat.GetComment() + "; fMet");
                        }
                        ChangeMade (CCleanupChange::eChangeComment);
                    }
                    return true;
                }
            }
        } else if (rna.IsSetExt() && !rna.GetExt().IsName()) {
            return false;
        } else if (type == CRNA_ref::eType_other) {
            // new convention follows ASN.1 spec comments, allows new RNA types
            return false;
        } else {
            // subsequent /product now added to comment
            string rna_name = "";
            string val = gb_qual.GetVal();
            if (rna.IsSetExt() && rna.GetExt().IsName()) {
                rna_name = rna.GetExt().GetName();
            }
            if (NStr::IsBlank (rna_name)) {
                rna.SetExt().SetName(val);
                ChangeMade (CCleanupChange::eChange_rRna);
                return true;
            } else {
                if (NStr::EqualNocase(rna_name, val)) {
                    return true;
                } else if (type == CRNA_ref::eType_miscRNA) {
                    return false;
                } else {
                    if (!feat.IsSetComment() || NStr::IsBlank (feat.GetComment())) {
                        feat.SetComment(val);
                    } else {
                        feat.SetComment (feat.GetComment() + "; " + val);
                    }
                    ChangeMade (CCleanupChange::eChangeComment);
                    return true;
                }
            }
        }

    }
    if (NStr::EqualNocase(qual, "anticodon")) {
        if (!rna.IsSetType()) {
            rna.SetType(CRNA_ref::eType_tRNA);
            ChangeMade(CCleanupChange::eChangeKeywords);
        }
        _ASSERT(rna.IsSetType());
        CRNA_ref::TType type = rna.GetType();
        if (type == CRNA_ref::eType_unknown) {
            rna.SetType(CRNA_ref::eType_tRNA);
            ChangeMade(CCleanupChange::eChangeKeywords);
        } else if (type != CRNA_ref::eType_tRNA) {
            return false;
        }
        if (!rna.IsSetExt()) {
            rna.SetExt().SetTRNA();
        }
        if ( rna.IsSetExt()  &&
             rna.GetExt().Which() == CRNA_ref::C_Ext::e_TRNA ) {
            
            CRef<CTrna_ext> trna = s_ParseTRnaFromAnticodonString (gb_qual.GetVal(), feat, m_Scope);
            if (trna->IsSetAa() || trna->IsSetAnticodon()) {
                /* don't apply at all if there are conflicts */
                bool apply_aa = false;
                bool apply_anticodon = false;
                bool ok_to_apply = true;
                
                /* look for conflict with aa */
                if (trna->IsSetAa()) {
                    if (rna.GetExt().GetTRNA().IsSetAa()) {
                        if (trna->GetAa().GetIupacaa() != rna.GetExt().GetTRNA().GetAa().GetIupacaa()) {
                            ok_to_apply = false;
                        }
                    } else {
                        apply_aa = true;
                    }
                }
                /* look for conflict with anticodon */
                if (trna->IsSetAnticodon()) {
                    if (rna.GetExt().GetTRNA().IsSetAnticodon()) {
                        if (sequence::Compare(rna.GetExt().GetTRNA().GetAnticodon(), trna->GetAnticodon(), m_Scope) != sequence::eSame) {
                            ok_to_apply = false;
                        }
                    } else {
                        apply_anticodon = true;
                    }
                }

                if (ok_to_apply) {
                    if (apply_aa) {
                        rna.SetExt().SetTRNA().SetAa().SetIupacaa(trna->GetAa().GetIupacaa());
                        ChangeMade (CCleanupChange::eChange_tRna);
                    }
                    if (apply_anticodon) {
                        CRef<CSeq_loc> anticodon(new CSeq_loc());
                        anticodon->Add (trna->GetAnticodon());
                        rna.SetExt().SetTRNA().SetAnticodon(*anticodon);
                        ChangeMade (CCleanupChange::eChangeAnticodon);
                    }
                    return true;
                }
            }
        }
    }
    return false;
}


bool CCleanup_imp::BasicCleanup(CProt_ref& prot, const CGb_qual& gb_qual)
{    
    const string& qual = gb_qual.GetQual();
    const string& val  = gb_qual.GetVal();

    if (NStr::EqualNocase(qual, "product")  ||  NStr::EqualNocase(qual, "standard_name")) {
        if (!prot.IsSetName()  ||  NStr::IsBlank(prot.GetName().front())) {
            prot.SetName().push_back(val);
            ChangeMade(CCleanupChange::eChangeQualifiers);
            if (prot.IsSetDesc()) {
                const CProt_ref::TDesc& desc = prot.GetDesc();
                FOR_EACH_NAME_ON_PROTREF (it, prot) {
                    if (NStr::EqualNocase(desc, *it)) {
                        prot.ResetDesc();
                        ChangeMade(CCleanupChange::eChangeQualifiers);
                        break;
                    }
                }
            }
            return true;
        }
    } else if (NStr::EqualNocase(qual, "function")) {
        prot.SetActivity().push_back(val);
        ChangeMade(CCleanupChange::eChangeQualifiers);
        return true;
    } else if (NStr::EqualNocase(qual, "EC_number")) {
        prot.SetEc().push_back(val);
        ChangeMade(CCleanupChange::eChangeQualifiers);
        return true;
    }

    return false;
}


static bool s_IsCompoundRptTypeValue( 
    const string& value )
//
//  Format of compound rpt_type values: (value[,value]*)
//
//  These are internal to sequin and are in theory cleaned up before the material
//  is released. However, some compound values have escaped into the wild and have 
//  not been retro-fixed yet (as of 2006-03-17).
//
{
    return ( NStr::StartsWith( value, '(' ) && NStr::EndsWith( value, ')' ) );
};


static void s_ExpandThisQual( 
    CSeq_feat::TQual& quals,        // the list of CGb_qual's.
    CSeq_feat::TQual::iterator& it, // points to the one qual we might expand.
    CSeq_feat::TQual& new_quals )    // new quals that will need to be inserted
//
//  Rules for "rpt_type" qualifiers (as of 2006-03-07):
//
//  There can be multiple occurrences of this qualifier, and we need to keep them 
//  all.
//  The value of this qualifier can also be a *list of values* which is *not* 
//  conforming to the ASN.1 and thus needs to be cleaned up. 
//
//  The cleanup entails turning the list of values into multiple occurrences of the 
//  given qualifier, each occurrence taking one of the values in the original 
//  list.
//
{
    CGb_qual& qual = **it;
    string  qual_type = qual.GetQual();
    string& val = qual.SetVal();
    if ( ! s_IsCompoundRptTypeValue( val ) ) {
        //
        //  nothing to do ...
        //
        return;
    }

    //
    //  Generate list of cleaned up values. Fix original qualifier and generate 
    //  list of new qualifiers to be added to the original list:
    //    
    vector< string > newValues;
    string valueList = val.substr(1, val.length() - 2);
    NStr::Tokenize(valueList, ",", newValues);
    
    qual.SetVal( newValues[0] );
    
    for ( size_t i=1; i < newValues.size(); ++i ) {
        CRef< CGb_qual > newQual( new CGb_qual() );
        newQual->SetQual( qual_type );
        newQual->SetVal( newValues[i] );
        new_quals.push_back( newQual ); 
    }
};

void CCleanup_imp::x_ExpandCombinedQuals(CSeq_feat::TQual& quals)
{
    CSeq_feat::TQual    new_quals;
    NON_CONST_ITERATE (CSeq_feat::TQual, it, quals) {
        CGb_qual& gb_qual = **it;
                
        string& qual = gb_qual.SetQual();
        
        if (NStr::EqualNocase(qual, "rpt_type")) {            
            s_ExpandThisQual( quals, it, new_quals ); 
        } else if (NStr::EqualNocase(qual, "rpt_unit")) {
            s_ExpandThisQual( quals, it, new_quals ); 
        } else if (NStr::EqualNocase(qual, "usedin")) {
            s_ExpandThisQual( quals, it, new_quals ); 
        } else if (NStr::EqualNocase(qual, "old_locus_tag")) {
            s_ExpandThisQual( quals, it, new_quals ); 
        } else if (NStr::EqualNocase(qual, "compare")) {
            s_ExpandThisQual( quals, it, new_quals ); 
        }
    }
    
    if ( ! new_quals.empty() ) {
        quals.insert(quals.end(), new_quals.begin(), new_quals.end());
        ChangeMade(CCleanupChange::eChangeQualifiers);
    }
}



bool CCleanup_imp::BasicCleanup(CSeq_feat& feat, CGb_qual& gb_qual)
{
    // return true if we should delete this gb_qual.
    
    _ASSERT(feat.IsSetQual()  &&  feat.IsSetData());

    CSeq_feat::TData& data = feat.SetData();
    string& qual = gb_qual.SetQual();
    string& val  = gb_qual.SetVal();

    if (NStr::EqualNocase (qual, "cons_splice")) {
        return true; // all cons_splice qualifiers should be removed
    }

    // 'replace' qualifier
    if (NStr::EqualNocase(qual, "replace")) {
        if (data.IsImp()  &&  data.GetImp().IsSetKey()) {
            CSeq_feat::TData::TImp& imp = feat.SetData().SetImp();
            if (NStr::EqualNocase(imp.GetKey(), "variation")) {
                NStr::ToLower(val);
            }
        }
        if ( ! NStr::IsBlank(val) &&  val.find_first_not_of("ACGTUacgtu") == NPOS) {
            NStr::ToLower(val);
            string val_no_u = NStr::Replace(val, "u", "t");
            if (val_no_u != val) {
                gb_qual.SetVal(val_no_u);
                ChangeMade(CCleanupChange::eChangeQualifiers);
            }
        }
    }
    
    if (NStr::EqualNocase(qual, "partial")) {
        feat.SetPartial();
        ChangeMade(CCleanupChange::eChangeQualifiers);
        return true;  // mark qual for deletion
    } else if (NStr::EqualNocase(qual, "evidence")) {
    /*
        if (NStr::EqualNocase(val, "experimental")) {
            if (!feat.IsSetExp_ev()  ||  feat.GetExp_ev() != CSeq_feat::eExp_ev_not_experimental) {
                feat.SetExp_ev(CSeq_feat::eExp_ev_experimental);
            }
        } else if (NStr::EqualNocase(val, "not_experimental")) {
            feat.SetExp_ev(CSeq_feat::eExp_ev_not_experimental);
        }
    */
        return true;  // mark qual for deletion
    } else if (NStr::EqualNocase(qual, "exception")) {
        feat.SetExcept(true);
        if (!NStr::IsBlank(val)  &&  !NStr::EqualNocase(val, "true")) {
            if (!feat.IsSetExcept_text()) {
                feat.SetExcept_text(val);
                ChangeMade(CCleanupChange::eChangeQualifiers);
            }
        }
        return true;  // mark qual for deletion
    } else if (NStr::EqualNocase(qual, "experiment")) {
        if (NStr::EqualNocase(val, "experimental evidence, no additional details recorded")) {
            ChangeMade(CCleanupChange::eChangeQualifiers);
            return true;  // mark qual for deletion
        }
    } else if (NStr::EqualNocase(qual, "inference")) {
        if (NStr::EqualNocase(val, "non-experimental evidence, no additional details recorded")) {
            ChangeMade(CCleanupChange::eChangeQualifiers);
            return true;  // mark qual for deletion
        }
    } else if (NStr::EqualNocase(qual, "note")  ||
               NStr::EqualNocase(qual, "notes")  ||
               NStr::EqualNocase(qual, "comment")) {
        if (!feat.IsSetComment()) {
            feat.SetComment(val);
        } else {
            (feat.SetComment() += "; ") += val;
        }
        ChangeMade(CCleanupChange::eChangeComment);
        ChangeMade(CCleanupChange::eChangeQualifiers);
        return true;  // mark qual for deletion
    } else if (NStr::EqualNocase(qual, "db_xref")) {
        string tag, db;
        if (NStr::SplitInTwo(val, ":", db, tag)) {
            CRef<CDbtag> dbp(new CDbtag);
            dbp->SetDb(db);
            dbp->SetTag().SetStr(tag);
            feat.SetDbxref().push_back(dbp);
            ChangeMade(CCleanupChange::eChangeDbxrefs);
            return true;  // mark qual for deletion
        }
    } else if (NStr::EqualNocase(qual, "gdb_xref")) {
        CRef<CDbtag> dbp(new CDbtag);
        dbp->SetDb("GDB");
        dbp->SetTag().SetStr(val);
        feat.SetDbxref().push_back(dbp);
        ChangeMade(CCleanupChange::eChangeDbxrefs);
        return true;  // mark qual for deletion
    } else if (NStr::EqualNocase(qual, "pseudo")) {
        feat.SetPseudo(true);
        ChangeMade(CCleanupChange::eChangeQualifiers);
        return true;  // mark qual for deletion
    } else if (data.IsGene()  &&  BasicCleanup(data.SetGene(), gb_qual)) {
        return true;  // mark qual for deletion
    } else if (data.IsCdregion()  &&  BasicCleanup(feat, data.SetCdregion(), gb_qual) ) {
        return true;  // mark qual for deletion
    } else if (data.IsRna()  &&  BasicCleanup(feat, data.SetRna(), gb_qual)) {
        return true;  // mark qual for deletion
    } else if (data.IsProt()  &&  BasicCleanup(data.SetProt(), gb_qual)) {
        return true;  // mark qual for deletion
    } else if (NStr::EqualNocase(qual, "gene")) {
        if (!NStr::IsBlank(val)) {
            CRef<CSeqFeatXref> xref(new CSeqFeatXref);
            xref->SetData().SetGene().SetLocus(val);
            feat.SetXref().push_back(xref);
            ChangeMade(CCleanupChange::eCopyGeneXref);
            return true;  // mark qual for deletion
        }
    } else if (NStr::EqualNocase(qual, "codon_start")) {
        if (!data.IsCdregion()) {
            // not legal on anything but CDS, so remove it
            return true;  // mark qual for deletion
        }
    }

    return false;
}


static bool s_IsJustQuotes(const string& str)
{
    ITERATE (string, it, str) {
        if ((*it > ' ')  &&  (*it != '"')  &&  (*it  != '\'')) {
            return false;
        }
    }
    return true;
}


void CCleanup_imp::BasicCleanup(CGb_qual& gbq)
{
    CLEAN_STRING_MEMBER_JUNK(gbq, Qual);
    if (!gbq.IsSetQual()) {
        gbq.SetQual(kEmptyStr);
    }
    
    CLEAN_STRING_MEMBER(gbq, Val);
    if (gbq.IsSetVal()  &&  s_IsJustQuotes(gbq.GetVal())) {
        gbq.ResetVal();
        ChangeMade(CCleanupChange::eCleanDoubleQuotes);
    }
    if (!gbq.IsSetVal()) {
        gbq.SetVal(kEmptyStr);
    }
    _ASSERT(gbq.IsSetQual()  &&  gbq.IsSetVal());
    
    if (NStr::EqualNocase(gbq.GetQual(), "cons_splice")) {
        x_CleanupConsSplice(gbq);
    } else if (NStr::EqualNocase(gbq.GetQual(), "rpt_unit")  ||
               NStr::EqualNocase(gbq.GetQual(), "rpt_unit_range")  ||
               NStr::EqualNocase(gbq.GetQual(), "rpt_unit_seq")) {
        bool range_qual = x_CleanupRptUnit(gbq);
        if (NStr::EqualNocase(gbq.GetQual(), "rpt_unit")) {
            if (range_qual) {
                gbq.SetQual("rpt_unit_range");
            } else {
                gbq.SetQual("rpt_unit_seq");
            }
            ChangeMade(CCleanupChange::eChangeQualifiers);
        }
    }
    x_ChangeTransposonToMobileElement(gbq);
    x_ChangeInsertionSeqToMobileElement(gbq);
}



// Gb_qual cleanup

void CCleanup_imp::x_ChangeInsertionSeqToMobileElement(CGb_qual& gbq)
//
//  As of Dec 2006, "insertion_seq" is no longer legal as a qualifier. The replacement
//  qualifier is "mobile_element". In addition, the value has to be massaged to
//  reflect the "insertion_seq".
//
{
    if (NStr::EqualNocase(gbq.GetQual(), "insertion_seq")) {
        gbq.SetQual("mobile_element");
        gbq.SetVal( string("insertion sequence: ") + gbq.GetVal() );
        ChangeMade(CCleanupChange::eChangeQualifiers);
    }
}


void CCleanup_imp::x_ChangeTransposonToMobileElement(CGb_qual& gbq)
//
//  As of Dec 2006, "transposon" is no longer legal as a qualifier. The replacement
//  qualifier is "mobile_element". In addition, the value has to be massaged to
//  indicate "integron" or "transposon".
//
{
    const string IntegronValues[] = {
        "class I integron",
        "class II integron",
        "class III integron",
        "class 1 integron",
        "class 2 integron",
        "class 3 integron"
    };
    const string* endIntegronValues 
        = IntegronValues + sizeof(IntegronValues)/sizeof(*IntegronValues);

    if (NStr::EqualNocase(gbq.GetQual(), "transposon")) {
        const string* pValue = std::find(IntegronValues, endIntegronValues, gbq.GetVal());

        gbq.SetQual("mobile_element");

        // If the value is one of the IntegronValues, change it to "integron: class XXX":

        if ( pValue != endIntegronValues ) {
            string::size_type cutoff = pValue->find( " integron" );
            _ASSERT( cutoff != string::npos ) /* typo in IntegronValues? */;
            gbq.SetVal( string("integron: ") + pValue->substr(0, cutoff) );
        }

        // Otherwise, just prefix it with "transposon: ":
        else {
            gbq.SetVal( string("transposon: ") + gbq.GetVal() );
        }
        
        ChangeMade(CCleanupChange::eChangeQualifiers);
    }
}


void CCleanup_imp::x_CleanupConsSplice(CGb_qual& gbq)

{
    string& val = gbq.SetVal();
    
    if (!NStr::StartsWith(val, "(5'site:")) {
        return;
    }
    
    size_t pos = val.find(",3'site:");
    if (pos != NPOS) {
        val.insert(pos + 1, " ");
        ChangeMade(CCleanupChange::eChangeQualifiers);
    }
}

static
bool s_HasUpper (string val)
{
    bool rval = false;

    string::iterator it = val.begin();
    while (it != val.end()) {
        if (isupper(*it)) {
            rval = true;
            break;
        }
        ++it;
    }
    return rval;
}


// return true if the val indicates this a range qualifier "[0-9]+..[0-9]+"

bool CCleanup_imp::x_CleanupRptUnit(CGb_qual& gbq)
{
    CGb_qual::TVal& val = gbq.SetVal();
    
    if (NStr::IsBlank(val)) {
        return false;
    }
    if( string::npos != val.find_first_not_of( "ACGTUNacgtun0123456789()" ) ) {
        if (s_HasUpper(val)) {
            val = NStr::ToLower(val);
            ChangeMade(CCleanupChange::eChangeQualifiers);
            return false;
        }
    } 
    bool    digits1, sep, digits2;
    digits1 = sep = digits2 = false;
    string s;
    string::const_iterator it = val.begin();
    string::const_iterator end = val.end();
    while (it != end) {
        while (it != end  &&  (*it == '('  ||  *it == ')'  ||  *it == ',')) {
            s += *it++;
        }
        while (it != end  &&  isspace((unsigned char)(*it))) {
            ++it;
        }
        while (it != end  &&  isdigit((unsigned char)(*it))) {
            s += *it++;
            digits1 = true;
        }
        if (it != end  &&  (*it == '.'  ||  *it == '-')) {
            while (it != end  &&  (*it == '.'  ||  *it == '-')) {
                ++it;
            }
            s += "..";
            sep = true;
        }
        while (it != end  &&  isspace((unsigned char)(*it))) {
            ++it;
        }
        while (it != end  &&  isdigit((unsigned char)(*it))) {
            s += *it++;
            digits2 = true;
        }
        while (it != end  &&  isspace((unsigned char)(*it))) {
            ++it;
        }
        if (it != end) {
            char c = *it;
            if (c != '('  &&  c != ')'  &&  c != ','  &&  c != '.'  &&
                !isspace((unsigned char) c)  &&  !isdigit((unsigned char) c)) {
                if (s_HasUpper(val)) {
                    val = NStr::ToLower(val);
                    ChangeMade(CCleanupChange::eChangeQualifiers);
                }
                return false;
            }
        }
    }
    if (val != s) {
        val = s;
        ChangeMade(CCleanupChange::eChangeQualifiers);
    }
    
    return  (digits1 && sep && digits2);
}


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE
