/* $Id: cleanup_general.cpp 161443 2009-05-27 20:45:44Z kans $
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
 *   Implementation of BasicCleanup for objects in objects/general
 *
 */

#include <ncbi_pch.hpp>
#include <objects/general/Object_id.hpp>
#include <objects/general/Dbtag.hpp>
#include <objects/general/Person_id.hpp>
#include <objects/general/Name_std.hpp>
#include <objects/general/User_object.hpp>
#include <objects/general/User_field.hpp>
#include <objects/general/Date_std.hpp>
#include <util/static_map.hpp>

#include <vector>

#include "cleanup_utils.hpp"
#include "cleanupp.hpp"


BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::



// convert str variant where the string is all digits to id variant
void CCleanup_imp::x_TagCleanup(CDbtag::TTag& tag)
{
    if (!tag.IsStr()) {
        return;
    }
    const CDbtag::TTag::TStr& str = tag.GetStr();
    if (NStr::IsBlank(str)) {
        return;
    }
    
    bool all_zero = true;
    ITERATE (CDbtag::TTag::TStr, it, str) {
        if (isdigit((unsigned char)(*it))) {
            if (*it != '0') {
                all_zero = false;
            }
        } else if (!isspace((unsigned char)(*it))) {
            return;
        }
    }
    
    if (str[0] != '0'  ||  all_zero) {
        try {
            tag.SetId(NStr::StringToUInt(str));
            ChangeMade(CCleanupChange::eChangeDbxrefs);
        } catch (CStringException&) {
            // just leave things as are
        }
    }
}


/** 
* replace old DB names with new ones:
*  Swiss-Prot       => UniProtKB/Swiss-Prot
*  SPTREMBL, TrEMBL => UniProtKB/TrEMBL
*  SUBTILIS         => SubtiList
*  LocusID          => GeneID
*  MaizeDB          => MaizeGDB
*  GeneW            => HGNC
*  MGD              => MGI
*  IFO              => NBRC
*  BHB              => BioHealthBase
**/
void CCleanup_imp::x_DbCleanup(string& db)
{
    if (NStr::EqualNocase(db, "Swiss-Prot")
        || NStr::EqualNocase (db, "SWISSPROT")) {
        db = "UniProtKB/Swiss-Prot";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "SPTREMBL")  ||
               NStr::EqualNocase(db, "TrEMBL")
               || NStr::EqualNocase(db, "UniProt/Swiss-prot")) {
        db = "UniProtKB/TrEMBL";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "SUBTILIS")) {
        db = "SubtiList";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "LocusID")) {
        db = "GeneID";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "MaizeDB")) {
        db = "MaizeGDB";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "GeneW")) {
        db = "HGNC";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "MGD")) {
        db = "MGI";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "IFO")) {
        db = "NBRC";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "BHB")) {
        db = "BioHealthBase";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase(db, "GeneDB")) {
        db = "BioHealthBase";
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    }
}


void CCleanup_imp::BasicCleanup(CDbtag& dbtag)
{
    if (dbtag.IsSetDb()) {
        if (CleanString(dbtag.SetDb(), true)) {
            ChangeMade(CCleanupChange::eTrimSpaces);
        }
    } else {
        dbtag.SetDb();  // make sure mandatory field is set.
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    }
    
    if (!dbtag.IsSetTag()) {
        dbtag.SetTag();  // make sure mandatory field is set.
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    }

    // need to do x_DbCleanup and x_TagCleanup so that
    // the steps afterward will be looking at the corrected db value
    // and the corrected tag value
    x_DbCleanup(dbtag.SetDb());
    x_TagCleanup(dbtag.SetTag());

    if (dbtag.GetTag().IsStr()) {
        if (CleanString(dbtag.SetTag().SetStr())) {
            ChangeMade(CCleanupChange::eTrimSpaces);
        }
    }
    if (NStr::EqualNocase(dbtag.GetDb(), "HPRD") && NStr::StartsWith (dbtag.GetTag().GetStr(), "HPRD_")) {
        dbtag.SetTag().SetStr (dbtag.GetTag().GetStr().substr (5));
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    } else if (NStr::EqualNocase (dbtag.GetDb(), "MGI") 
               && (NStr::StartsWith (dbtag.GetTag().GetStr(), "MGI:")
                   || NStr::StartsWith (dbtag.GetTag().GetStr(), "MGD:"))) {
        dbtag.SetTag().SetStr (dbtag.GetTag().GetStr().substr (4));
        ChangeMade(CCleanupChange::eChangeDbxrefs);
    }

    _ASSERT(dbtag.IsSetDb()  &&  dbtag.IsSetTag());  // mandatory fields
    
    dbtag.InvalidateType();
}


void CCleanup_imp::BasicCleanup(CPerson_id& pid, bool fix_initials)
{
    switch (pid.Which()) {
        case CPerson_id::e_Name:
            BasicCleanup(pid.SetName(), fix_initials);
            break;
        case CPerson_id::e_Ml:
            TRUNCATE_CHOICE_SPACES(pid, Ml);
            break;
        case CPerson_id::e_Str:
            TRUNCATE_CHOICE_SPACES(pid, Str);
            break;
        case CPerson_id::e_Consortium:
            TRUNCATE_CHOICE_SPACES(pid, Consortium);
            break;
        default:
            break;
    }
}



/*****  CName_std   Basic Cleanup  */


void CCleanup_imp::BasicCleanup(CName_std& name, bool fix_initials)
{
    // before cleanup, get information from initials
    if (name.IsSetInitials()) {
        if (!name.IsSetSuffix()) {
            x_ExtractSuffixFromInitials(name);
        }
        x_FixEtAl(name);
    }

    // do strings cleanup
    CLEAN_STRING_MEMBER(name, Last);
    CLEAN_STRING_MEMBER(name, First);
    CLEAN_STRING_MEMBER(name, Middle);
    CLEAN_STRING_MEMBER(name, Full);
    CLEAN_STRING_MEMBER(name, Initials);
    CLEAN_STRING_MEMBER(name, Suffix);
    CLEAN_STRING_MEMBER(name, Title);

    if (name.IsSetSuffix()) {
        x_FixSuffix(name);
    }

    // regenerate initials from first name
    if (fix_initials) {
        x_FixInitials(name);
    }
}


void CCleanup_imp::x_FixEtAl(CName_std& name)
{
    _ASSERT(name.IsSetInitials());

    if (name.IsSetLast()  &&  NStr::Equal(name.GetLast(), "et")  &&
        (NStr::Equal(name.GetInitials(), "al")  ||  NStr::Equal(name.GetInitials(), "al."))  &&
        (!name.IsSetFirst()  ||  NStr::IsBlank(name.GetFirst()))) {
        name.ResetFirst();
        name.ResetInitials();
        name.SetLast("et al.");
        ChangeMade(CCleanupChange::eNormalizeAuthors);
    }
}


static bool s_FixInitials(string& initials)
{
    string fixed;

    string::iterator it = initials.begin();
    while (it != initials.end()  &&  !isalpha((unsigned char)(*it))) {  // skip junk
        ++it;
    }
    char prev = '\0';
    while (it != initials.end()) {
        char c = *it;
    
        if (c == ',') {
            c = '.';
        }
        if (isalpha((unsigned char) c)) {
            if (prev != '\0'  &&  isupper((unsigned char) c)  &&  isupper((unsigned char) prev)) {
                fixed += '.';
            }
        } else if (c == '-') {
            if (prev != '\0'  &&  prev != '.') {
                fixed += '.';
            }
        } else if (c == '.') {
            if (prev == '-'  ||  prev == '.') {
                ++it;
                continue;
            }
        } else {
            ++it;
            continue;
        }
        fixed += c;

        prev = c;
        ++it;
    }
    if (!fixed.empty()  &&  fixed[fixed.length() - 1] == '-') {
        fixed.resize(fixed.length() - 1);
    }
    if (initials != fixed) {
        initials = fixed;
        return true;
    }
    return false;
}


CName_std::TInitials s_GetInitialsFromFirst(const CName_std& name)
{
    _ASSERT(name.IsSetFirst());
    
    const CName_std::TFirst& first = name.GetFirst();
    if (NStr::IsBlank(first)) {
        return kEmptyStr;
    }
    
    bool up = name.IsSetInitials()  &&  
        !name.GetInitials().empty()  &&  
        isupper((unsigned char) name.GetInitials()[0]);
    
    string initials;
    string::const_iterator end = first.end();
    string::const_iterator it = first.begin();
    while (it != end) {
        // skip "junk"
        while (it != end  &&  !isalpha((unsigned char)(*it))) {
            switch( (unsigned char)(*it) ) {

            case '(':
                // regard words in parentheses as nick names. Nick names do not contribute
                // to the initials
                //
                while ( ')' != (unsigned char)(*it)  &&  it != end ) {
                    ++it;
                }
                ++it;
                break;

            default:
                ++it;
                break;

            }
        }
        // copy the first character of each word
        if (it != end) {
            if (!initials.empty()) {
                char c = *initials.rbegin();
                if (isupper((unsigned char)(*it))  &&  c != '.'  &&  c != '-') {
                    initials += '.';
                }
            }
            initials += up ? toupper((unsigned char)(*it)) : *it;
        }
        // skip the rest of the word
        while (it != end  &&  !isspace((unsigned char)(*it))  &&  *it != ','  &&  *it != '-') {
            ++it;
        }
        if (it != end  &&  *it == '-') {
            initials += '.';
            initials += *it++;
        }
        up = false;
    }
    if (!initials.empty()  &&  !(*initials.rbegin() == '-')) {
        initials += '.';
    }
    
    return initials;
}


void CCleanup_imp::x_FixInitials(CName_std& name)
{
    if (name.IsSetInitials()) {
        if (s_FixInitials(name.SetInitials())) {
            ChangeMade(CCleanupChange::eNormalizeAuthors);
        }
        
        if (name.SetInitials().empty()) {
            name.ResetInitials();
            ChangeMade(CCleanupChange::eNormalizeAuthors);
        }
    }

    if (name.IsSetFirst()  &&  !name.IsSetInitials()) {
        string new_initials = s_GetInitialsFromFirst(name);
        if ( ! new_initials.empty() ) {
            name.SetInitials(new_initials);
            ChangeMade(CCleanupChange::eNormalizeAuthors);
        }
        return;
    }
    if (!name.IsSetInitials()) {
        return;
    }

    if (! name.IsSetSuffix()) {
        CName_std::TInitials& initials = name.SetInitials();
        if (NStr::EndsWith (initials, ".Jr.")) {
            initials.erase (initials.end() - 3);
            name.SetSuffix ("Jr.");
        } else if (NStr::EndsWith (initials, ".Sr.")) {
            initials.erase (initials.end() - 3);
            name.SetSuffix ("Sr.");
        }
    }

    if (!name.IsSetFirst()) {
        return;
    }
    _ASSERT(name.IsSetFirst()  &&  name.IsSetInitials());

    CName_std::TInitials& initials = name.SetInitials();
    string f_initials = s_GetInitialsFromFirst(name);
    if (initials == f_initials) {
        return;
    }
    
    string::iterator it1 = f_initials.begin();
    string::iterator it2 = initials.begin();
    while (it1 != f_initials.end()  &&  it2 != initials.end()  &&
        toupper((unsigned char)(*it1)) == toupper((unsigned char)(*it2))) {
        ++it1;
        ++it2;
    }
    if (it2 != initials.end()  &&  it1 != f_initials.end()  &&  *it1 == '.'  &&  islower((unsigned char)(*it2))) {
        f_initials.erase(it1);
    }
    while (it2 != initials.end()) {
        f_initials += *it2++;
    }
    if (f_initials != initials) {
        initials.swap(f_initials);
        ChangeMade(CCleanupChange::eNormalizeAuthors);
    }
}


// mapping of wrong suffixes to the correct ones.
typedef pair<string, string> TStringPair;
static const TStringPair bad_sfxs[] = {
    TStringPair("1d"  , "I"),
    TStringPair("1st" , "I"),
    TStringPair("2d"  , "II"),
    TStringPair("2nd" , "II"),
    TStringPair("3d"  , "III"),
    TStringPair("3rd" , "III"),
    TStringPair("4th" , "IV"),
    TStringPair("5th" , "V"),
    TStringPair("6th" , "VI"),
    //TStringPair("I."  , "I"),
    TStringPair("II." , "II"),
    TStringPair("III.", "III"),
    TStringPair("IV." , "IV"),
    TStringPair("Jr"  , "Jr."),
    TStringPair("Sr"  , "Sr."),    
    //TStringPair("V."  , "V"),
    TStringPair("VI." , "VI")
};
typedef CStaticArrayMap<string, string> TSuffixMap;
DEFINE_STATIC_ARRAY_MAP(TSuffixMap, sc_BadSuffixes, bad_sfxs);

// move the suffix from the initials field to the suffix field.
void CCleanup_imp::x_ExtractSuffixFromInitials(CName_std& name)
{
    _ASSERT(name.IsSetInitials()  &&  !name.IsSetSuffix());

    string& initials = name.SetInitials();
    size_t len = initials.length();

    if (initials.find('.') == NPOS) {
        return;
    }

    // look for standard suffixes in intials
    typedef CName_std::TSuffixes TSuffixes;
    const TSuffixes& suffixes = CName_std::GetStandardSuffixes();
    TSuffixes::const_iterator best = suffixes.end();
    ITERATE (TSuffixes, it, suffixes) {
        if (NStr::EndsWith(initials, *it)) {
            if (best == suffixes.end()) {
                best = it;
            } else if (best->length() < it->length()) {
                best = it;
            }
        }
    }
    if (best != suffixes.end()) {
        initials.resize(len - best->length());
        name.SetSuffix(*best);
        ChangeMade(CCleanupChange::eNormalizeAuthors);
        return;
    }

    // look for bad suffixes in initials
    ITERATE (TSuffixMap, it, sc_BadSuffixes) {
        if (NStr::EndsWith(initials, it->first)  &&
            initials.length() > it->first.length()) {
            initials.resize(len - it->first.length());
            name.SetSuffix(it->second);
            ChangeMade(CCleanupChange::eNormalizeAuthors);
            return;
        }
    }
}


void CCleanup_imp::x_FixSuffix(CName_std& name)
{
    _ASSERT(name.IsSetSuffix());

    CName_std::TSuffix& suffix = name.SetSuffix();

    ITERATE (TSuffixMap, it, sc_BadSuffixes) {
        if (NStr::EqualNocase(suffix, it->first)) {
            name.SetSuffix(it->second);
            ChangeMade(CCleanupChange::eNormalizeAuthors);
            break;
        }
    }
}



// === User Objects & Fields

void CCleanup_imp::x_CleanupUserString(string& str)
{
    if (!NStr::IsBlank(str)) {
        size_t str_n = str.size();
        NStr::TruncateSpacesInPlace(str);
        if (str_n != str.size()) {
            ChangeMade(CCleanupChange::eTrimSpaces);
        }
    }
}

void CCleanup_imp::BasicCleanup(CObject_id& oid)
{
    if (oid.IsStr()) {
        x_CleanupUserString(oid.SetStr());
    }
}


void CCleanup_imp::BasicCleanup(CUser_field& field)
{
    if (field.IsSetLabel()) {
        BasicCleanup(field.SetLabel());
    }
    
    if (field.IsSetData()) {
        CUser_field::TData& data = field.SetData();
        switch (data.Which()) {
            case CUser_field::TData::e_Str:
                x_CleanupUserString(data.SetStr());
                break;
            case CUser_field::TData::e_Strs:
                NON_CONST_ITERATE (CUser_field::TData::TStrs, it, data.SetStrs()) {
                    x_CleanupUserString(*it);
                }
                break;
            case CUser_field::TData::e_Object:
                BasicCleanup(data.SetObject());
                break;
            case CUser_field::TData::e_Objects:
                NON_CONST_ITERATE (CUser_field::TData::TObjects, it, data.SetObjects()) {
                    BasicCleanup(**it);
                }
                break;
            case CUser_field::TData::e_Fields:
                NON_CONST_ITERATE (CUser_field::TData::TFields, it, data.SetFields()) {
                    BasicCleanup(**it);
                }
                break;
            default:
                break;
        };
    }
}


void CCleanup_imp::BasicCleanup(CUser_object& uo)
{
    if (uo.IsSetType()) {
        BasicCleanup(uo.SetType());
    }
    
    NON_CONST_ITERATE (CUser_object::TData, it, uo.SetData()) {
        BasicCleanup(**it);
    }
}


void CCleanup_imp::BasicCleanup(CDate_std& date)
{
    if (date.IsSetSecond() && date.GetSecond() > 59) {
        date.ResetSecond();
        ChangeMade(CCleanupChange::eCleanupDate);
    }

    if (date.IsSetMinute() && date.GetMinute() > 59) {
        date.ResetMinute();
        date.ResetSecond();
        ChangeMade(CCleanupChange::eCleanupDate);
    }
    
    if (date.IsSetHour() && date.GetHour() > 23) {
        date.ResetHour();
        date.ResetMinute();
        date.ResetSecond();
        ChangeMade(CCleanupChange::eCleanupDate);
    }
}


void CCleanup_imp::BasicCleanup(CDate& date)
{
    if (date.IsStd()) {
        BasicCleanup(date.SetStd());
    }
}


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE
