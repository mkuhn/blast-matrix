/* $Id: RNA_ref.cpp 170000 2009-09-08 13:00:17Z bollin $
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
#include <objects/seqfeat/RNA_ref.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

// destructor
CRNA_ref::~CRNA_ref(void)
{
}


typedef pair<const string, CRNA_ref::TType> TRnaTypePair;
static const TRnaTypePair sc_rna_type_map[] = {
    TRnaTypePair("mRNA", CRNA_ref::eType_mRNA),
    TRnaTypePair("misc_RNA", CRNA_ref::eType_other),
    TRnaTypePair("ncRNA", CRNA_ref::eType_ncRNA),
    TRnaTypePair("precursor_RNA", CRNA_ref::eType_premsg),
    TRnaTypePair("rRNA", CRNA_ref::eType_rRNA),
    TRnaTypePair("scRNA", CRNA_ref::eType_scRNA),
    TRnaTypePair("snRNA", CRNA_ref::eType_snRNA),
    TRnaTypePair("snoRNA", CRNA_ref::eType_snoRNA),
    TRnaTypePair("tRNA", CRNA_ref::eType_tRNA),
    TRnaTypePair("tmRNA", CRNA_ref::eType_tmRNA)
};

typedef CStaticArrayMap<const string, CRNA_ref::TType> TRnaTypeMap;
DEFINE_STATIC_ARRAY_MAP(TRnaTypeMap, sc_RnaTypeMap, sc_rna_type_map);

string CRNA_ref::GetRnaTypeName (const CRNA_ref::EType rna_type)
{
    string rna_type_name = "";
    TRnaTypeMap::const_iterator rna_type_it = sc_RnaTypeMap.begin();
    while (rna_type_it != sc_RnaTypeMap.end() && rna_type_it->second != rna_type) {
        ++rna_type_it;
    }
    if (rna_type_it != sc_RnaTypeMap.end()) {
        rna_type_name = rna_type_it->first;
    }
    return rna_type_name;
}



END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE

/* Original file checksum: lines: 65, chars: 1885, CRC32: 6e6c3f8a */
