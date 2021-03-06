/* $Id: Dbtag.cpp 240007 2011-02-02 16:24:38Z rafanovi $
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
 *   using specifications from the ASN data definition file
 *   'general.asn'.
 *
 * ---------------------------------------------------------------------------
 */

// standard includes

// generated includes
#include <ncbi_pch.hpp>
#include <objects/general/Dbtag.hpp>
#include <objects/general/Object_id.hpp>
#include <corelib/ncbistd.hpp>
#include <util/static_map.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

// When adding to these lists, please take care to keep them in
// case-sensitive sorted order (lowercase entries last).

typedef pair<const char*, CDbtag::EDbtagType> TDbxrefPair;
static const TDbxrefPair kApprovedDbXrefs[] = {
    TDbxrefPair("AFTOL", CDbtag::eDbtagType_AFTOL),
    TDbxrefPair("APHIDBASE", CDbtag::eDbtagType_APHIDBASE),
    TDbxrefPair("ASAP", CDbtag::eDbtagType_ASAP),
    TDbxrefPair("ATCC", CDbtag::eDbtagType_ATCC),
    TDbxrefPair("ATCC(dna)", CDbtag::eDbtagType_ATCC_dna),
    TDbxrefPair("ATCC(in host)", CDbtag::eDbtagType_ATCC_in_host),
    TDbxrefPair("AceView/WormGenes", CDbtag::eDbtagType_AceView_WormGenes),
    TDbxrefPair("AntWeb", CDbtag::eDbtagType_AntWeb),
    TDbxrefPair("ApiDB", CDbtag::eDbtagType_ApiDB),
    TDbxrefPair("ApiDB_CryptoDB", CDbtag::eDbtagType_ApiDB_CryptoDB),
    TDbxrefPair("ApiDB_PlasmoDB", CDbtag::eDbtagType_ApiDB_PlasmoDB),
    TDbxrefPair("ApiDB_ToxoDB", CDbtag::eDbtagType_ApiDB_ToxoDB),
    TDbxrefPair("Axeldb", CDbtag::eDbtagType_axeldb),
    TDbxrefPair("BDGP_EST", CDbtag::eDbtagType_BDGP_EST),
    TDbxrefPair("BDGP_INS", CDbtag::eDbtagType_BDGP_INS),
    TDbxrefPair("BEETLEBASE", CDbtag::eDbtagType_BEETLEBASE),
    TDbxrefPair("BOLD", CDbtag::eDbtagType_BoLD),
    TDbxrefPair("CDD", CDbtag::eDbtagType_CDD),
    TDbxrefPair("CK", CDbtag::eDbtagType_CK),
    TDbxrefPair("COG", CDbtag::eDbtagType_COG),
    TDbxrefPair("ENSEMBL", CDbtag::eDbtagType_ENSEMBL),
    TDbxrefPair("ERIC", CDbtag::eDbtagType_ERIC),
    TDbxrefPair("ESTLIB", CDbtag::eDbtagType_ESTLIB),
    TDbxrefPair("EcoGene", CDbtag::eDbtagType_EcoGene),
    TDbxrefPair("FANTOM_DB", CDbtag::eDbtagType_FANTOM_DB),
    TDbxrefPair("FLYBASE", CDbtag::eDbtagType_FLYBASE),
    TDbxrefPair("GABI", CDbtag::eDbtagType_GABI),
    TDbxrefPair("GDB", CDbtag::eDbtagType_GDB),
    TDbxrefPair("GI", CDbtag::eDbtagType_GI),
    TDbxrefPair("GO", CDbtag::eDbtagType_GO),
    TDbxrefPair("GOA", CDbtag::eDbtagType_GOA),
    TDbxrefPair("GRIN", CDbtag::eDbtagType_GRIN),
    TDbxrefPair("GeneDB", CDbtag::eDbtagType_GeneDB),
    TDbxrefPair("GeneID", CDbtag::eDbtagType_GeneID),
    TDbxrefPair("Greengenes", CDbtag::eDbtagType_Greengenes),
    TDbxrefPair("H-InvDB", CDbtag::eDbtagType_H_InvDB),
    TDbxrefPair("HGNC", CDbtag::eDbtagType_HGNC),
    TDbxrefPair("HMP", CDbtag::eDbtagType_HMP),
    TDbxrefPair("HOMD", CDbtag::eDbtagType_HOMD),
    TDbxrefPair("HSSP", CDbtag::eDbtagType_HSSP),
    TDbxrefPair("IMGT/GENE-DB", CDbtag::eDbtagType_IMGT_GENEDB),
    TDbxrefPair("IMGT/HLA", CDbtag::eDbtagType_IMGT_HLA),
    TDbxrefPair("IMGT/LIGM", CDbtag::eDbtagType_IMGT_LIGM),
    TDbxrefPair("IRD", CDbtag::eDbtagType_IRD),
    TDbxrefPair("ISD", CDbtag::eDbtagType_ISD),
    TDbxrefPair("ISFinder", CDbtag::eDbtagType_ISFinder),
    TDbxrefPair("InterPro", CDbtag::eDbtagType_Interpro),
    TDbxrefPair("InterimID", CDbtag::eDbtagType_InterimID),
    TDbxrefPair("JCM", CDbtag::eDbtagType_JCM),
    TDbxrefPair("JGIDB", CDbtag::eDbtagType_JGIDB),
    TDbxrefPair("LocusID", CDbtag::eDbtagType_LocusID),
    TDbxrefPair("MGI", CDbtag::eDbtagType_MGI),
    TDbxrefPair("MIM", CDbtag::eDbtagType_MIM),
    TDbxrefPair("MaizeGDB", CDbtag::eDbtagType_MaizeGDB),
    TDbxrefPair("MycoBank", CDbtag::eDbtagType_MycoBank),
    TDbxrefPair("NBRC", CDbtag::eDbtagType_IFO),
    TDbxrefPair("NMPDR", CDbtag::eDbtagType_NMPDR),
    TDbxrefPair("NRESTdb", CDbtag::eDbtagType_NRESTdb),
    TDbxrefPair("NextDB", CDbtag::eDbtagType_NextDB),
    TDbxrefPair("Osa1", CDbtag::eDbtagType_Osa1),
    TDbxrefPair("PBmice", CDbtag::eDbtagType_PBmice),
    TDbxrefPair("PDB", CDbtag::eDbtagType_PDB),
    TDbxrefPair("PFAM", CDbtag::eDbtagType_PFAM),
    TDbxrefPair("PGN", CDbtag::eDbtagType_PGN),
    TDbxrefPair("PIR", CDbtag::eDbtagType_PIR),
    TDbxrefPair("PSEUDO", CDbtag::eDbtagType_PSEUDO),
    TDbxrefPair("Pathema", CDbtag::eDbtagType_Pathema),
    TDbxrefPair("PseudoCap", CDbtag::eDbtagType_PseudoCap),
    TDbxrefPair("RAP-DB", CDbtag::eDbtagType_RAP_DB),
    TDbxrefPair("RATMAP", CDbtag::eDbtagType_RATMAP),
    TDbxrefPair("RFAM", CDbtag::eDbtagType_RFAM),
    TDbxrefPair("RGD", CDbtag::eDbtagType_RGD),
    TDbxrefPair("RZPD", CDbtag::eDbtagType_RZPD),
    TDbxrefPair("RiceGenes", CDbtag::eDbtagType_RiceGenes),
    TDbxrefPair("SEED", CDbtag::eDbtagType_SEED),
    TDbxrefPair("SGD", CDbtag::eDbtagType_SGD),
    TDbxrefPair("SGN", CDbtag::eDbtagType_SGN),
    TDbxrefPair("SoyBase", CDbtag::eDbtagType_SoyBase),
    TDbxrefPair("SubtiList", CDbtag::eDbtagType_SubtiList),
    TDbxrefPair("TIGRFAM", CDbtag::eDbtagType_TIGRFAM),
    TDbxrefPair("UNILIB", CDbtag::eDbtagType_UNILIB),
    TDbxrefPair("UNITE", CDbtag::eDbtagType_UNITE),
    TDbxrefPair("UniGene", CDbtag::eDbtagType_UniGene),
    TDbxrefPair("UniProtKB/Swiss-Prot", CDbtag::eDbtagType_UniProt_SwissProt),
    TDbxrefPair("UniProtKB/TrEMBL", CDbtag::eDbtagType_UniProt_TrEMBL),
    TDbxrefPair("UniSTS", CDbtag::eDbtagType_UniSTS),
    TDbxrefPair("VBASE2", CDbtag::eDbtagType_VBASE2),
    TDbxrefPair("VectorBase", CDbtag::eDbtagType_VectorBase),
    TDbxrefPair("WorfDB", CDbtag::eDbtagType_WorfDB),
    TDbxrefPair("WormBase", CDbtag::eDbtagType_WormBase),
    TDbxrefPair("Xenbase", CDbtag::eDbtagType_Xenbase),
    TDbxrefPair("ZFIN", CDbtag::eDbtagType_ZFIN),
    TDbxrefPair("dbClone", CDbtag::eDbtagType_dbClone),
    TDbxrefPair("dbCloneLib", CDbtag::eDbtagType_dbCloneLib),
    TDbxrefPair("dbEST", CDbtag::eDbtagType_dbEST),
    TDbxrefPair("dbProbe", CDbtag::eDbtagType_dbProbe),
    TDbxrefPair("dbSNP", CDbtag::eDbtagType_dbSNP),
    TDbxrefPair("dbSTS", CDbtag::eDbtagType_dbSTS),
    TDbxrefPair("dictyBase", CDbtag::eDbtagType_dictyBase),
    TDbxrefPair("niaEST", CDbtag::eDbtagType_niaEST), 
    TDbxrefPair("taxon", CDbtag::eDbtagType_taxon)
};

static const TDbxrefPair kApprovedRefSeqDbXrefs[] = {
    TDbxrefPair("BEEBASE", CDbtag::eDbtagType_BEEBASE),
    TDbxrefPair("BioProject", CDbtag::eDbtagType_BioProject),
    TDbxrefPair("CCDS", CDbtag::eDbtagType_CCDS),
    TDbxrefPair("CGNC", CDbtag::eDbtagType_CGNC),
    TDbxrefPair("CloneID", CDbtag::eDbtagType_CloneID),
    TDbxrefPair("DDBJ", CDbtag::eDbtagType_DDBJ),
    TDbxrefPair("ECOCYC", CDbtag::eDbtagType_ECOCYC),
    TDbxrefPair("EMBL", CDbtag::eDbtagType_EMBL),
    TDbxrefPair("HPRD", CDbtag::eDbtagType_HPRD),
    TDbxrefPair("LRG", CDbtag::eDbtagType_LRG),
    TDbxrefPair("NASONIABASE", CDbtag::eDbtagType_NASONIABASE),
    TDbxrefPair("PBR", CDbtag::eDbtagType_PBR),
    TDbxrefPair("REBASE", CDbtag::eDbtagType_REBASE),
    TDbxrefPair("SK-FST", CDbtag::eDbtagType_SK_FST),
    TDbxrefPair("TAIR", CDbtag::eDbtagType_TAIR),
    TDbxrefPair("VBRC", CDbtag::eDbtagType_VBRC),
    TDbxrefPair("miRBase", CDbtag::eDbtagType_miRBase)
};

static const TDbxrefPair kApprovedSrcDbXrefs[] = {
    TDbxrefPair("AFTOL", CDbtag::eDbtagType_AFTOL),
    TDbxrefPair("ATCC", CDbtag::eDbtagType_ATCC),
    TDbxrefPair("ATCC(dna)", CDbtag::eDbtagType_ATCC_dna),
    TDbxrefPair("ATCC(in host)", CDbtag::eDbtagType_ATCC_in_host),
    TDbxrefPair("AntWeb", CDbtag::eDbtagType_AntWeb),
    TDbxrefPair("BOLD", CDbtag::eDbtagType_BoLD),
    TDbxrefPair("FANTOM_DB", CDbtag::eDbtagType_FANTOM_DB),
    TDbxrefPair("FLYBASE", CDbtag::eDbtagType_FLYBASE),
    TDbxrefPair("GRIN", CDbtag::eDbtagType_GRIN),
    TDbxrefPair("HMP", CDbtag::eDbtagType_HMP),
    TDbxrefPair("HOMD", CDbtag::eDbtagType_HOMD),
    TDbxrefPair("IMGT/HLA", CDbtag::eDbtagType_IMGT_HLA),
    TDbxrefPair("IMGT/LIGM", CDbtag::eDbtagType_IMGT_LIGM),
    TDbxrefPair("JCM", CDbtag::eDbtagType_JCM),
    TDbxrefPair("MGI", CDbtag::eDbtagType_MGI),
    TDbxrefPair("MycoBank", CDbtag::eDbtagType_MycoBank),
    TDbxrefPair("NBRC", CDbtag::eDbtagType_IFO),
    TDbxrefPair("RZPD", CDbtag::eDbtagType_RZPD),
    TDbxrefPair("UNILIB", CDbtag::eDbtagType_UNILIB),
    TDbxrefPair("UNITE", CDbtag::eDbtagType_UNITE), 
    TDbxrefPair("taxon", CDbtag::eDbtagType_taxon)
};

static const TDbxrefPair kApprovedProbeDbXrefs[] = {
    TDbxrefPair("GEO", CDbtag::eDbtagType_GEO)
};


static const char* const kSkippableDbXrefs[] = {
    "BankIt",
    "NCBIFILE",
    "TMSMART"
};

// case sensetive
typedef CStaticArrayMap<const char*, CDbtag::EDbtagType, PCase_CStr> TDbxrefTypeMap;
// case insensitive, per the C Toolkit
typedef CStaticArraySet<const char*, PNocase_CStr> TDbxrefSet;

DEFINE_STATIC_ARRAY_MAP(TDbxrefTypeMap, sc_ApprovedDb,       kApprovedDbXrefs);
DEFINE_STATIC_ARRAY_MAP(TDbxrefTypeMap, sc_ApprovedRefSeqDb, kApprovedRefSeqDbXrefs);
DEFINE_STATIC_ARRAY_MAP(TDbxrefTypeMap, sc_ApprovedSrcDb,    kApprovedSrcDbXrefs);
DEFINE_STATIC_ARRAY_MAP(TDbxrefTypeMap, sc_ApprovedProbeDb,  kApprovedProbeDbXrefs);
DEFINE_STATIC_ARRAY_MAP(TDbxrefSet,     sc_SkippableDbXrefs, kSkippableDbXrefs);

// destructor
CDbtag::~CDbtag(void)
{
}

bool CDbtag::Match(const CDbtag& dbt2) const
{
    if (! PNocase().Equals(GetDb(), dbt2.GetDb()))
        return false;
    return ((GetTag()).Match((dbt2.GetTag())));
}


int CDbtag::Compare(const CDbtag& dbt2) const
{
    int ret = PNocase().Compare(GetDb(), dbt2.GetDb());
    if (ret == 0) {
        ret = GetTag().Compare(dbt2.GetTag());
    }
    return ret;
}


// Appends a label to "label" based on content of CDbtag 
void CDbtag::GetLabel(string* label) const
{
    const CObject_id& id = GetTag();
    switch (id.Which()) {
    case CObject_id::e_Str:
        *label += GetDb() + ": " + id.GetStr();
        break;
    case CObject_id::e_Id:
        *label += GetDb() + ": " + NStr::IntToString(id.GetId());
        break;
    default:
        *label += GetDb();
    }
}

// Test if CDbtag.db is in the approved databases list.
// NOTE: 'GenBank', 'EMBL', 'DDBJ' and 'REBASE' are approved only in 
//        the context of a RefSeq record.
bool CDbtag::IsApproved( EIsRefseq refseq, EIsSource is_source, EIsEstOrGss is_est_or_gss ) const
{
    if ( !CanGetDb() ) {
        return false;
    }
    const string& db = GetDb();

    if( refseq == eIsRefseq_Yes && sc_ApprovedRefSeqDb.find(db.c_str()) != sc_ApprovedRefSeqDb.end() ) {
        return true;
    }

    if( is_source == eIsSource_Yes ) {
        bool found = ( sc_ApprovedSrcDb.find(db.c_str()) != sc_ApprovedSrcDb.end() );
        if ( ! found && (is_est_or_gss == eIsEstOrGss_Yes) ) {
            // special case: for EST or GSS, source features are allowed non-src dbxrefs
            found = ( sc_ApprovedDb.find(db.c_str()) != sc_ApprovedDb.end() );
        }
        return found;
    } else {
        return sc_ApprovedDb.find(db.c_str()) != sc_ApprovedDb.end();
    }
}


const char* CDbtag::IsApprovedNoCase(EIsRefseq refseq, EIsSource is_source ) const
{
    if ( !CanGetDb() ) {
        return NULL;
    }
    const string& db = GetDb();
    
    const char* retval = NULL;
    // This is *slow*.  Someone needs to replace this with a binary search or something
    // Since this function isn't even used right now, I'm postponing fixing this.
    ITERATE (TDbxrefTypeMap, it, sc_ApprovedDb) {
        if ( NStr::EqualNocase(db, it->first) ) {
            retval = it->first;
            break;
        }
    }
    // This is *slow*.  Someone needs to replace this with a binary search or something
    // Since this function isn't even used right now, I'm postponing fixing this.
    if ( retval == NULL  &&  (refseq == eIsRefseq_Yes) ) {
        ITERATE (TDbxrefTypeMap, it, sc_ApprovedRefSeqDb) {
            if ( NStr::EqualNocase(db, it->first) ) {
                retval = it->first;
                break;
            }
        }
    }

    return retval;
}


bool CDbtag::IsApproved(TDbtagGroup group) const
{
    if ( !CanGetDb() ) {
        return false;
    }
    const string& db = GetDb();

    if ( (group & fGenBank) != 0 && sc_ApprovedDb.find(db.c_str()) != sc_ApprovedDb.end()) {
        return true;
    }

    if ( (group & fRefSeq) != 0 && sc_ApprovedRefSeqDb.find(db.c_str()) != sc_ApprovedRefSeqDb.end()) {
        return true;
    }

    if ( (group & fSrc) != 0 && sc_ApprovedSrcDb.find(db.c_str()) != sc_ApprovedSrcDb.end()) {
        return true;
    }

    if ( (group & fProbe) != 0 && sc_ApprovedProbeDb.find(db.c_str()) != sc_ApprovedProbeDb.end()) {
        return true;
    }

    return false;
}


bool CDbtag::IsSkippable(void) const
{
    return sc_SkippableDbXrefs.find(GetDb().c_str())
        != sc_SkippableDbXrefs.end();
}


// Retrieve the enumerated type for the dbtag
CDbtag::EDbtagType CDbtag::GetType(void) const
{
    if (m_Type == eDbtagType_bad) {
        if ( !CanGetDb() ) {
            return m_Type;
        }

        const string& db = GetDb();

        TDbxrefTypeMap::const_iterator iter;

        iter = sc_ApprovedDb.find(db.c_str());
        if ( iter != sc_ApprovedDb.end() ) {
            m_Type = iter->second;
            return m_Type;
        }

        iter = sc_ApprovedRefSeqDb.find(db.c_str());
        if ( iter != sc_ApprovedRefSeqDb.end() ) {
            m_Type = iter->second;
            return m_Type;
        }

        iter = sc_ApprovedSrcDb.find(db.c_str());
        if ( iter != sc_ApprovedSrcDb.end() ) {
            m_Type = iter->second;
            return m_Type;
        }

        iter = sc_ApprovedProbeDb.find(db.c_str());
        if ( iter != sc_ApprovedProbeDb.end() ) {
            m_Type = iter->second;
            return m_Type;
        }
    }

    return m_Type;
}


CDbtag::TDbtagGroup CDbtag::GetDBFlags (string& correct_caps) const
{
    correct_caps = "";
    CDbtag::TDbtagGroup rsult = fNone;

    if ( !CanGetDb() ) {
        return fNone;
    }
    const string& db = GetDb();
    
    ITERATE (TDbxrefTypeMap, it, sc_ApprovedDb) {
        if ( NStr::EqualNocase(db, it->first) ) {
            correct_caps = it->first;
            rsult |= fGenBank;
        }
    }

    ITERATE (TDbxrefTypeMap, it, sc_ApprovedRefSeqDb) {
        if ( NStr::EqualNocase(db, it->first) ) {
            correct_caps = it->first;
            rsult |= fRefSeq;
        }
    }

    ITERATE (TDbxrefTypeMap, it, sc_ApprovedSrcDb) {
        if ( NStr::EqualNocase(db, it->first) ) {
            correct_caps = it->first;
            rsult |= fSrc;
        }
    }

    ITERATE (TDbxrefTypeMap, it, sc_ApprovedProbeDb) {
        if ( NStr::EqualNocase(db, it->first) ) {
            correct_caps = it->first;
            rsult |= fProbe;
        }
    }

    return rsult;
}


bool CDbtag::GetDBFlags (bool& is_refseq, bool& is_src, string& correct_caps) const
{
    CDbtag::TDbtagGroup group = CDbtag::GetDBFlags(correct_caps);

    is_refseq = ((group & fRefSeq) != 0);
    is_src    = ((group & fSrc)    != 0);

    return group != fNone;
}


// Force a refresh of the internal type
void CDbtag::InvalidateType(void)
{
    m_Type = eDbtagType_bad;
}


//=========================================================================//
//                              URLs                                       //
//=========================================================================//

// special case URLs
static const string kFBan = "http://www.fruitfly.org/cgi-bin/annot/fban?";
static const string kHInvDbHIT = "http://www.jbirc.aist.go.jp/hinv/hinvsys/servlet/ExecServlet?KEN_INDEX=0&KEN_TYPE=30&KEN_STR=";
static const string kHInvDbHIX = "http://www.jbirc.aist.go.jp/hinv/hinvsys/servlet/ExecServlet?KEN_INDEX=0&KEN_TYPE=31&KEN_STR=";
static const string kDictyPrim = "http://dictybase.org/db/cgi-bin/gene_page.pl?primary_id=";
static const string kMiRBaseMat = "http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=";

// mapping of DB to its URL; please sort these by tag name (mostly in
// case-sensitive ASCII-betical order as above)
typedef pair<CDbtag::EDbtagType, string>    TDbtUrl;
static const TDbtUrl sc_url_prefix[] = {
    TDbtUrl(CDbtag::eDbtagType_AFTOL, "http://aftol1.biology.duke.edu/pub/displayTaxonInfo?aftol_id="),
    TDbtUrl(CDbtag::eDbtagType_APHIDBASE, "http://webapps1.genouest.org/grs-1.0/grs?reportID=chado_genome_report&objectID="),
    TDbtUrl(CDbtag::eDbtagType_ASAP, "https://asap.ahabs.wisc.edu/annotation/php/feature_info.php?FeatureID="),
    TDbtUrl(CDbtag::eDbtagType_ATCC, "http://www.atcc.org/SearchCatalogs/linkin?id="),
    TDbtUrl(CDbtag::eDbtagType_AceView_WormGenes, "http://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/av.cgi?db=worm&c=gene&q="),
    TDbtUrl(CDbtag::eDbtagType_AntWeb, "http://www.antweb.org/specimen.do?name="),
    TDbtUrl(CDbtag::eDbtagType_ApiDB, "http://www.apidb.org/apidb/showRecord.do?name=GeneRecordClasses.ApiDBGeneRecordClass&primary_key="),
    TDbtUrl(CDbtag::eDbtagType_ApiDB_CryptoDB, "http://cryptodb.org/cryptodb/showRecord.do?name=GeneRecordClasses.GeneRecordClass&project_id=&primary_key="),
    TDbtUrl(CDbtag::eDbtagType_ApiDB_PlasmoDB, "http://www.plasmodb.org/plasmo/showRecord.do?name=GeneRecordClasses.GeneRecordClass&project_id=&primary_key="),
    TDbtUrl(CDbtag::eDbtagType_ApiDB_ToxoDB, "http://www.toxodb.org/toxo/showRecord.do?name=GeneRecordClasses.GeneRecordClass&project_id=&primary_key="),
    TDbtUrl(CDbtag::eDbtagType_BEETLEBASE, "http://www.beetlebase.org/cgi-bin/report.cgi?name="),
    TDbtUrl(CDbtag::eDbtagType_BoLD, "http://www.boldsystems.org/connectivity/specimenlookup.php?processid="),
    TDbtUrl(CDbtag::eDbtagType_CCDS, "http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA="),
    TDbtUrl(CDbtag::eDbtagType_CDD, "http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="),
    TDbtUrl(CDbtag::eDbtagType_CK, "http://flybane.berkeley.edu/cgi-bin/cDNA/CK_clone.pl?db=CK&dbid="),
    TDbtUrl(CDbtag::eDbtagType_COG, "http://www.ncbi.nlm.nih.gov/COG/new/release/cow.cgi?cog="),
    TDbtUrl(CDbtag::eDbtagType_ECOCYC, "http://biocyc.org/ECOLI/new-image?type=GENE&object="),
    TDbtUrl(CDbtag::eDbtagType_ENSEMBL, "http://www.ensembl.org/id/"),
    TDbtUrl(CDbtag::eDbtagType_ERIC, "http://www.ericbrc.org/genbank/dbxref/"),
    TDbtUrl(CDbtag::eDbtagType_EcoGene, "http://ecogene.org/geneInfo.php?eg_id="),
    TDbtUrl(CDbtag::eDbtagType_FANTOM_DB, "http://fantom.gsc.riken.jp/db/annotate/main.cgi?masterid="),
    TDbtUrl(CDbtag::eDbtagType_FLYBASE, "http://flybase.bio.indiana.edu/.bin/fbidq.html?"),
    TDbtUrl(CDbtag::eDbtagType_GABI, "http://www.gabipd.org/database/cgi-bin/GreenCards.pl.cgi?Mode=ShowSequence&App=ncbi&SequenceId="),
    TDbtUrl(CDbtag::eDbtagType_GO, "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&depth=1&query=GO:"),
    TDbtUrl(CDbtag::eDbtagType_GOA, "http://www.ebi.ac.uk/ego/GProtein?ac="),
    TDbtUrl(CDbtag::eDbtagType_GRIN, "http://www.ars-grin.gov/cgi-bin/npgs/acc/display.pl?"),
    TDbtUrl(CDbtag::eDbtagType_GeneDB, "http://old.genedb.org/genedb/Search?organism=All%3A*&name="),
    TDbtUrl(CDbtag::eDbtagType_GeneID, "http://www.ncbi.nlm.nih.gov/gene/"),
    TDbtUrl(CDbtag::eDbtagType_Greengenes, "http://greengenes.lbl.gov/cgi-bin/show_one_record_v2.pl?prokMSA_id="),
    TDbtUrl(CDbtag::eDbtagType_HGNC, "http://www.genenames.org/data/hgnc_data.php?hgnc_id="),
    TDbtUrl(CDbtag::eDbtagType_HMP, "http://www.hmpdacc-resources.org/cgi-bin/hmp_catalog/main.cgi?section=HmpSummary&page=displayHmpProject&hmp_id="),
    TDbtUrl(CDbtag::eDbtagType_HOMD, "http://www.homd.org/"),
    TDbtUrl(CDbtag::eDbtagType_HPRD, "http://www.hprd.org/protein/"),
    TDbtUrl(CDbtag::eDbtagType_HSSP, "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-newId+-e+hssp-ID:"),
    TDbtUrl(CDbtag::eDbtagType_H_InvDB, "http://www.h-invitational.jp"),
    TDbtUrl(CDbtag::eDbtagType_IFO, "http://www.nbrc.nite.go.jp/NBRC2/NBRCCatalogueDetailServlet?ID=NBRC&CAT="),
    TDbtUrl(CDbtag::eDbtagType_IMGT_GENEDB, "http://imgt.cines.fr/cgi-bin/GENElect.jv?species=Homo+sapiens&query=2+"),
    TDbtUrl(CDbtag::eDbtagType_IMGT_LIGM, "http://imgt.cines.fr:8104/cgi-bin/IMGTlect.jv?query=202+"),
    TDbtUrl(CDbtag::eDbtagType_IRD, "http://www.fludb.org/brc/fluSegmentDetails.do?irdSubmissionId="),
    TDbtUrl(CDbtag::eDbtagType_ISD, "http://www.flu.lanl.gov/search/view_record.html?accession="),
    TDbtUrl(CDbtag::eDbtagType_ISFinder, "http://www-is.biotoul.fr/scripts/is/is_spec.idc?name="),
    TDbtUrl(CDbtag::eDbtagType_InterimID, "http://www.ncbi.nlm.nih.gov/gene/"),
    TDbtUrl(CDbtag::eDbtagType_Interpro, "http://www.ebi.ac.uk/interpro/ISearch?mode=ipr&query="),
    TDbtUrl(CDbtag::eDbtagType_JCM, "http://www.jcm.riken.go.jp/cgi-bin/jcm/jcm_number?JCM="),
    TDbtUrl(CDbtag::eDbtagType_JGIDB, "http://genome.jgi-psf.org/cgi-bin/jgrs?id="),
    TDbtUrl(CDbtag::eDbtagType_LocusID, "http://www.ncbi.nlm.nih.gov/gene/"),
    TDbtUrl(CDbtag::eDbtagType_MGI, "http://www.informatics.jax.org/searches/accession_report.cgi?id=MGI:"),
    TDbtUrl(CDbtag::eDbtagType_MIM, "http://www.ncbi.nlm.nih.gov/omim/"),
    TDbtUrl(CDbtag::eDbtagType_MaizeGDB, "http://www.maizegdb.org/cgi-bin/displaylocusrecord.cgi?id="),
    TDbtUrl(CDbtag::eDbtagType_MycoBank, "http://www.mycobank.org/MycoTaxo.aspx?Link=T&Rec="),
    TDbtUrl(CDbtag::eDbtagType_NMPDR, "http://www.nmpdr.org/linkin.cgi?id="),
    TDbtUrl(CDbtag::eDbtagType_NRESTdb, "http://genome.ukm.my/nrestdb/db/single_view_est.php?id="),
    TDbtUrl(CDbtag::eDbtagType_NextDB, "http://nematode.lab.nig.ac.jp/cgi-bin/db/ShowGeneInfo.sh?celk="),
    TDbtUrl(CDbtag::eDbtagType_Osa1, "http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice/?name="),
    TDbtUrl(CDbtag::eDbtagType_PBR, "http://www.poxvirus.org/query.asp?web_id="),
    TDbtUrl(CDbtag::eDbtagType_PBmice, "http://www.idmshanghai.cn/PBmice/DetailedSearch.do?type=insert&id="),
    TDbtUrl(CDbtag::eDbtagType_PDB, "http://www.rcsb.org/pdb/cgi/explore.cgi?pdbId="),
    TDbtUrl(CDbtag::eDbtagType_PFAM, "http://pfam.sanger.ac.uk/family?acc="),
    TDbtUrl(CDbtag::eDbtagType_PGN, "http://pgn.cornell.edu/cgi-bin/search/seq_search_result.pl?identifier="),
    TDbtUrl(CDbtag::eDbtagType_Pathema, "http://pathema.jcvi.org/cgi-bin/Burkholderia/shared/GenePage.cgi?all=1&locus="),
    TDbtUrl(CDbtag::eDbtagType_PseudoCap, "http://www.pseudomonas.com/getAnnotation.do?locusID="),
    TDbtUrl(CDbtag::eDbtagType_RAP_DB, "http://rapdb.dna.affrc.go.jp/cgi-bin/gbrowse_details/latest?name="),
    TDbtUrl(CDbtag::eDbtagType_RATMAP, "http://ratmap.gen.gu.se/ShowSingleLocus.htm?accno="),
    TDbtUrl(CDbtag::eDbtagType_REBASE, "http://rebase.neb.com/rebase/enz/"),
    TDbtUrl(CDbtag::eDbtagType_RFAM, "http://www.sanger.ac.uk/cgi-bin/Rfam/getacc?"),
    TDbtUrl(CDbtag::eDbtagType_RGD, "http://rgd.mcw.edu/query/query.cgi?id="),
    TDbtUrl(CDbtag::eDbtagType_RiceGenes, "http://ars-genome.cornell.edu/cgi-bin/WebAce/webace?db=ricegenes&class=Marker&object="),
    TDbtUrl(CDbtag::eDbtagType_SEED, "http://www.theseed.org/linkin.cgi?id="),
    TDbtUrl(CDbtag::eDbtagType_SGD, "http://db.yeastgenome.org/cgi-bin/SGD/locus.pl?locus="),
    TDbtUrl(CDbtag::eDbtagType_SGN, "http://www.sgn.cornell.edu/search/est.pl?request_type=7&request_id="),
    TDbtUrl(CDbtag::eDbtagType_SK_FST, "http://aafc-aac.usask.ca/fst/"),
    TDbtUrl(CDbtag::eDbtagType_SubtiList, "http://genolist.pasteur.fr/SubtiList/genome.cgi?external_query+"),
    TDbtUrl(CDbtag::eDbtagType_TAIR, "http://www.arabidopsis.org/servlets/TairObject?type=locus&name="),
    TDbtUrl(CDbtag::eDbtagType_TIGRFAM, "http://cmr.tigr.org/tigr-scripts/CMR/HmmReport.cgi?hmm_acc="),
    TDbtUrl(CDbtag::eDbtagType_UNITE, "http://unite.ut.ee/bl_forw.php?nimi="),
    TDbtUrl(CDbtag::eDbtagType_UniGene, "http://www.ncbi.nlm.nih.gov/unigene?term="),
    TDbtUrl(CDbtag::eDbtagType_UniProt_SwissProt, "http://www.uniprot.org/uniprot/"),
    TDbtUrl(CDbtag::eDbtagType_UniProt_TrEMBL, "http://www.uniprot.org/uniprot/"),
    TDbtUrl(CDbtag::eDbtagType_UniSTS, "http://www.ncbi.nlm.nih.gov/genome/sts/sts.cgi?uid="),
    TDbtUrl(CDbtag::eDbtagType_VBASE2, "http://www.dnaplot.de/vbase2/vgene.php?id="),
    TDbtUrl(CDbtag::eDbtagType_VBRC, "http://vbrc.org/query.asp?web_view=curation&web_id="),
    TDbtUrl(CDbtag::eDbtagType_VectorBase, "http://www.vectorbase.org/Genome/BRCGene/?feature="),
    TDbtUrl(CDbtag::eDbtagType_WorfDB, "http://worfdb.dfci.harvard.edu/search.pl?form=1&search="),
    TDbtUrl(CDbtag::eDbtagType_WormBase, "http://www.wormbase.org/db/gene/gene?class=CDS;name="),
    TDbtUrl(CDbtag::eDbtagType_Xenbase, "http://www.xenbase.org/gene/showgene.do?method=display&geneId="),
    TDbtUrl(CDbtag::eDbtagType_ZFIN, "http://zfin.org/cgi-bin/webdriver?MIval=aa-markerview.apg&OID="),
    TDbtUrl(CDbtag::eDbtagType_axeldb, "http://www.dkfz-heidelberg.de/tbi/services/axeldb/clone/xenopus?name="),
    TDbtUrl(CDbtag::eDbtagType_dbClone, "http://www.ncbi.nlm.nih.gov/sites/entrez?db=clone&cmd=Retrieve&list_uids="),
    TDbtUrl(CDbtag::eDbtagType_dbCloneLib, "http://www.ncbi.nlm.nih.gov/sites/entrez?db=clonelib&cmd=Retrieve&list_uids="),
    TDbtUrl(CDbtag::eDbtagType_dbEST, "http://www.ncbi.nlm.nih.gov/nucest/"),
    TDbtUrl(CDbtag::eDbtagType_dbProbe, "http://www.ncbi.nlm.nih.gov/sites/entrez?db=probe&cmd=Retrieve&list_uids="),
    TDbtUrl(CDbtag::eDbtagType_dbSNP, "http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs="),
    TDbtUrl(CDbtag::eDbtagType_dbSTS, "http://www.ncbi.nlm.nih.gov/nuccore/"),
    TDbtUrl(CDbtag::eDbtagType_dictyBase, "http://dictybase.org/db/cgi-bin/gene_page.pl?dictybaseid="),
    TDbtUrl(CDbtag::eDbtagType_miRBase, "http://microrna.sanger.ac.uk/cgi-bin/sequences/mirna_entry.pl?acc="),
    TDbtUrl(CDbtag::eDbtagType_niaEST, "http://lgsun.grc.nia.nih.gov/cgi-bin/pro3?sname1="),
    TDbtUrl(CDbtag::eDbtagType_taxon, "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?"),
    TDbtUrl(CDbtag::eDbtagType_BEEBASE, "http://genomes.arc.georgetown.edu/cgi-bin/gbrowse/bee_genome4/?name="),
    TDbtUrl(CDbtag::eDbtagType_NASONIABASE, "http://genomes.arc.georgetown.edu/cgi-bin/gbrowse/nasonia10_scaffold/?name=")
};

typedef CStaticArrayMap<CDbtag::EDbtagType, string> TUrlPrefixMap;
DEFINE_STATIC_ARRAY_MAP(TUrlPrefixMap, sc_UrlMap, sc_url_prefix);


string CDbtag::GetUrl(void) const
{
    TUrlPrefixMap::const_iterator it = sc_UrlMap.find(GetType());
    if (it == sc_UrlMap.end()) {
        return kEmptyStr;
    }
    const string* prefix = &(it->second);

    string tag;
    if (GetTag().IsStr()) {
        tag = GetTag().GetStr();
    } else if (GetTag().IsId()) {
        tag = NStr::IntToString(GetTag().GetId());
    }
    if (NStr::IsBlank(tag)) {
        return kEmptyStr;
    }

    // URLs are constructed by catenating the URL prefix with the specific tag
    // except in a few cases handled below.
    switch (GetType()) {
        case CDbtag::eDbtagType_FLYBASE:
            if (NStr::Find(tag, "FBan") != NPOS) {
                prefix = &kFBan;
            }
            break;

        case eDbtagType_MGI:
        case eDbtagType_MGD:
            if (NStr::StartsWith(tag, "MGI:", NStr::eNocase)  ||
                NStr::StartsWith(tag, "MGD:", NStr::eNocase)) {
                tag = tag.substr(4);
            }
            break;

        case eDbtagType_PID:
            if (tag[0] == 'g') {
                tag = tag.substr(1);
            }
            break;

        case eDbtagType_dbEST:
            tag.insert(0, "val=gnl|dbest|");
            break;

        case eDbtagType_dbSTS:
            break;

        case eDbtagType_niaEST:
            tag += "&val=1";
            break;

        case eDbtagType_MaizeGDB:
        case eDbtagType_IFO:
            tag.erase();
            break;

        case eDbtagType_GDB:
        {{
            SIZE_TYPE pos = NStr::Find(tag, "G00-");
            if (pos != NPOS) {
                tag = tag.substr(pos + 4);
                remove(tag.begin(), tag.end(), '-');
            } else if (!isdigit((unsigned char) tag[0])) {
                return kEmptyStr;
            }
            break;
        }}
        case eDbtagType_REBASE:
            tag += ".html";
            break;
        case eDbtagType_H_InvDB:
            if (NStr::Find(tag, "HIT")) {
                prefix = &kHInvDbHIT;
            } else if (NStr::Find(tag, "HIX")) {
                prefix = &kHInvDbHIX;
            }
            break;

        case eDbtagType_SK_FST:
            return *prefix;
            break;

        case CDbtag::eDbtagType_taxon:
            if (isdigit((unsigned char) tag[0])) {
                tag.insert(0, "id=");
            } else {
                tag.insert(0, "name=");
            }
            break;

        case CDbtag::eDbtagType_dictyBase:
            if (NStr::Find(tag, "_") != NPOS) {
                prefix = &kDictyPrim;
            }
            break;


        case CDbtag::eDbtagType_miRBase:
            if (NStr::Find(tag, "MIMAT") != NPOS) {
                prefix = &kMiRBaseMat;
            }
            break;

        case CDbtag::eDbtagType_GeneDB:
            tag +="&organism=pombe";
            break;

        default:
            break;
    }

    return *prefix + tag;
}


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE
