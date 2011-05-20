/* $Id: Cit_gen.cpp 254755 2011-02-16 21:42:42Z dicuccio $
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
 *   using specifications from the data definition file
 *   'biblio.asn'.
 */

// standard includes

// generated includes
#include <ncbi_pch.hpp>
#include <objects/biblio/Cit_gen.hpp>
#include <objects/biblio/label_util.hpp>
#include <objects/general/Date.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

// destructor
CCit_gen::~CCit_gen(void)
{
}


void CCit_gen::GetLabel(string* label, bool unique) const
{
    if (IsSetSerial_number()) {
        *label += "[" + NStr::IntToString(GetSerial_number()) + "]";
    }
    if (IsSetMuid()) {
        *label += "NLM" + NStr::IntToString(GetMuid());
    }

    string date;
    string* date_ptr = 0;
    if ( IsSetDate() ) {
        date_ptr = &date;
        GetDate().GetDate(date_ptr, true);
    }

    const string* title2 = 0;
    const string* titleunique = 0;
    bool unpublished = false;
    const CTitle* title = IsSetJournal() ? &GetJournal() : 0;
    const CAuth_list* authors = IsSetAuthors() ? &GetAuthors() : 0;
    const string* volume = IsSetVolume() ? &GetVolume() : 0;
    const string* issue = IsSetIssue() ? &GetIssue() : 0;
    const string* pages = IsSetPages() ? &GetPages() : 0;

    if (IsSetCit()) {
        if ( NStr::EqualNocase( GetCit(), "Unpublished") ) {
            unpublished = true;
        } else if (!title) {
            title2 = &GetCit();
        }
    }
    if (IsSetTitle()) {
        titleunique = &GetTitle();
    } else if (title2) {
        titleunique = title2;
    } else if (!title && IsSetCit()) {
        titleunique = &GetCit();
    }
    if (!title && !authors && !IsSetTitle() && !volume &&
        !pages && !issue) {
        titleunique = 0;
        if (IsSetCit()) {
            string cit(GetCit());
            if (!unique) {
                try {
                    cit.resize(cit.find_last_of('|'));
                } catch(length_error&) {}
            }   
            *label += cit;
        }
        return;
    }

    GetLabelContent(label, unique,
        authors, 0, title, 0, 0, 0, title2, titleunique,
        date_ptr, volume, issue, pages, unpublished);
}


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE

/* Original file checksum: lines: 64, chars: 1875, CRC32: 5ca91cd9 */
